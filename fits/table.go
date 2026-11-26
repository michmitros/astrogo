package fits

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Table is a FITS table extension (ASCII TABLE or BINTABLE) parsed into typed
// columns. Use ToTable to create a Table from a decoded HDU.
type Table struct {
	// Type is the extension type: ExtTable or ExtBinTable.
	Type ExtensionType
	// ExtName is the EXTNAME keyword value, empty if the keyword is absent.
	ExtName string
	// RowCount is the number of rows (NAXIS2).
	RowCount int
	// Columns holds the columns in the order they appear in the header.
	Columns []Column
	// colIndex maps column names (TTYPEn) to indexes in Columns.
	// If a name appears more than once, the last column wins.
	colIndex map[string]int
	// Headers contains the raw header cards of the table HDU.
	Headers Headers
}

// ColumnMeta holds the metadata of a single table column, read from the TTYPEn,
// TFORMn, TUNITn, TDIMn, TSCALn, TZEROn, TNULLn and TBCOLn keywords. Scale and
// Zero are exposed for the caller to apply; column data is returned unscaled.
type ColumnMeta struct {
	Index    int     // 0-based column index
	Name     string  // TTYPEn
	Format   string  // TFORMn
	Unit     string  // TUNITn
	Dim      []int   // TDIMn
	Scale    float64 // TSCALn (default 1.0)
	Zero     float64 // TZEROn (default 0.0)
	NullVal  string  // TNULLn
	AsciiCol int     // TBCOLn (1-based, ASCII tables only; 0 for BINTABLE)
}

// Column is a single column of a Table. The concrete type is always a
// *TypedColumn[T], where T depends on the column format and repeat count.
type Column interface {
	// Meta returns the column metadata.
	Meta() ColumnMeta
	// Len returns the number of rows in the column.
	Len() int
	// GetAny returns the value at the given 0-based row, or nil if the value
	// is undefined (null). It panics if row is out of range.
	GetAny(row int) any
}

// TypedColumn is a column of values of type T, one value per row. Null marks
// undefined values: when Null is non-nil and Null[row] is true, the value at
// that row is undefined and GetAny returns nil. Null is nil for columns that
// cannot hold undefined values (strings, bit arrays and multi-value columns).
type TypedColumn[T any] struct {
	meta ColumnMeta
	Data []T
	Null []bool
}

// Meta returns the column metadata.
func (c *TypedColumn[T]) Meta() ColumnMeta {
	return c.meta
}

// Len returns the number of rows in the column.
func (c *TypedColumn[T]) Len() int {
	return len(c.Data)
}

// GetAny returns the value at the given 0-based row, or nil if the value is
// undefined (null). It panics if row is out of range.
func (c *TypedColumn[T]) GetAny(row int) any {
	if c.Null != nil && c.Null[row] {
		return nil
	}
	return c.Data[row]
}

// BitArray holds the value of a BINTABLE bit column ('X' format): BitCount
// bits packed into Bytes, most significant bit of the first byte first.
// Bits beyond BitCount in the last byte are padding and carry no meaning.
type BitArray struct {
	BitCount int
	Bytes    []byte
}

// GetColumn returns the data of the named column (TTYPEn) as a []T. T must
// match the column's Go type exactly, e.g. int32 for format "1J", []int32 for
// "10J", string for "A" formats and BitArray for "X". It returns an error if
// the column does not exist or has a different type. The returned slice is
// the column's backing data, not a copy.
func GetColumn[T any](t *Table, name string) ([]T, error) {
	idx, ok := t.colIndex[name]
	if !ok {
		return nil, fmt.Errorf("column %s not found", name)
	}

	col, ok := t.Columns[idx].(*TypedColumn[T])
	if !ok {
		return nil, fmt.Errorf("invalid type for column %s", name)
	}

	return col.Data, nil
}

// ToTable parses a TABLE (ASCII) or BINTABLE extension HDU into a Table.
//
// ASCII columns map to string ('A'), int64 ('I') and float64 ('F', 'E', 'D').
// BINTABLE columns map to bool ('L'), uint8 ('B'), int16 ('I'), int32 ('J'),
// int64 ('K'), float32 ('E'), float64 ('D'), complex64 ('C'), complex128 ('M'),
// string ('A') and BitArray ('X'). Columns with a repeat count greater than 1
// (except 'A' and 'X') hold []T values per row. Variable-length array columns
// ('P', 'Q') hold []T values read from the heap that follows the row data in
// HDU.Data; only the numeric and logical element types are supported.
//
// Undefined values are reported through TypedColumn.Null: blank ASCII fields,
// ASCII integers equal to TNULLn, binary integers equal to TNULLn, NaN floats
// and logical bytes other than 'T' and 'F'. TSCALn and TZEROn are stored in
// ColumnMeta but not applied to the data.
//
// An error is returned if the HDU is not a table extension, if required
// keywords are missing or invalid, or if the data is too short for the
// declared rows and heap.
func ToTable(u HDU) (*Table, error) {
	switch u.Extension() {
	case ExtTable:
		return toTable(u)
	case ExtBinTable:
		return toBinTable(u)
	default:
		return nil, fmt.Errorf(
			"hdu is not a table: expected TABLE or BINTABLE, got XTENSION=%s", u.Extension(),
		)
	}
}

func toTable(u HDU) (*Table, error) {
	rowWidth, rowCount, nfields, err := tableParams(u)
	if err != nil {
		return nil, err
	}

	if pc := parseIntHeader(u.Headers["PCOUNT"]); pc != 0 {
		return nil, fmt.Errorf("invalid PCOUNT for TABLE: %d, expected 0", pc)
	}

	t := newTable(u, rowCount, nfields)
	for n := 1; n <= nfields; n++ {
		meta, err := columnMeta(u.Headers, n)
		if err != nil {
			return nil, err
		}

		code, width, err := parseAsciiFormat(meta.Format)
		if err != nil {
			return nil, fmt.Errorf("column %d: %w", n, err)
		}

		bcol, err := intHeader(u.Headers, fmt.Sprintf("TBCOL%d", n))
		if err != nil {
			return nil, err
		}

		if bcol < 1 || bcol-1+width > rowWidth {
			return nil, fmt.Errorf(
				"column %d: TBCOL%d=%d with width %d exceeds row width %d",
				n, n, bcol, width, rowWidth,
			)
		}

		meta.AsciiCol = bcol

		col, err := asciiColumn(u.Data, rowWidth, rowCount, bcol-1, width, code, meta)
		if err != nil {
			return nil, fmt.Errorf("column %d (%s): %w", n, meta.Name, err)
		}

		t.addColumn(col)
	}

	return t, nil
}

func toBinTable(u HDU) (*Table, error) {
	rowWidth, rowCount, nfields, err := tableParams(u)
	if err != nil {
		return nil, err
	}

	pcount := parseIntHeader(u.Headers["PCOUNT"])
	if pcount < 0 {
		return nil, fmt.Errorf("invalid PCOUNT: %d", pcount)
	}

	heapStart := rowWidth * rowCount

	if len(u.Data) < heapStart+pcount {
		return nil, fmt.Errorf(
			"data too short: have %d bytes, need %d (including %d heap bytes)",
			len(u.Data), heapStart+pcount, pcount,
		)
	}

	t := newTable(u, rowCount, nfields)

	off := 0
	for n := 1; n <= nfields; n++ {
		meta, err := columnMeta(u.Headers, n)
		if err != nil {
			return nil, err
		}

		repeat, code, elem, err := parseBinFormat(meta.Format)
		if err != nil {
			return nil, fmt.Errorf("column %d: %w", n, err)
		}

		width, err := binColumnWidth(repeat, code)
		if err != nil {
			return nil, fmt.Errorf("column %d: %w", n, err)
		}

		if off+width > rowWidth {
			return nil, fmt.Errorf(
				"column %d: width %d at offset %d exceeds row width %d",
				n, width, off, rowWidth,
			)
		}

		col, err := binColumn(u, meta, rowWidth, rowCount, off, repeat, code, elem, heapStart)
		if err != nil {
			return nil, fmt.Errorf("column %d (%s): %w", n, meta.Name, err)
		}

		t.addColumn(col)
		off += width
	}

	return t, nil
}

func tableParams(u HDU) (rowWidth, rowCount, nfields int, err error) {
	if naxis := u.Naxis(); naxis != 2 {
		return 0, 0, 0, fmt.Errorf("expected NAXIS=2 for a table, got %d", naxis)
	}

	if bitpix := u.Bitpix(); bitpix != 8 {
		return 0, 0, 0, fmt.Errorf("expected BITPIX=8 for a table, got %d", bitpix)
	}

	rowWidth, err = intHeader(u.Headers, "NAXIS1")
	if err != nil {
		return
	}

	if rowWidth <= 0 {
		err = fmt.Errorf("invalid NAXIS1: %d", rowWidth)
		return
	}

	rowCount, err = intHeader(u.Headers, "NAXIS2")
	if err != nil {
		return
	}

	if rowCount < 0 {
		err = fmt.Errorf("invalid NAXIS2: %d", rowCount)
		return
	}

	nfields, err = intHeader(u.Headers, "TFIELDS")
	if err != nil {
		return
	}

	if nfields < 0 {
		err = fmt.Errorf("invalid TFIELDS: %d", nfields)
		return
	}

	if len(u.Data) < rowWidth*rowCount {
		err = fmt.Errorf("data too short: have %d bytes, need %d", len(u.Data), rowWidth*rowCount)
		return
	}

	return
}

func newTable(u HDU, rowCount, nfields int) *Table {
	return &Table{
		Type:     u.Extension(),
		ExtName:  u.ExtName(),
		RowCount: rowCount,
		Columns:  make([]Column, 0, nfields),
		colIndex: make(map[string]int, nfields),
		Headers:  u.Headers,
	}
}

func (t *Table) addColumn(c Column) {
	idx := len(t.Columns)
	t.Columns = append(t.Columns, c)

	if name := c.Meta().Name; name != "" {
		t.colIndex[name] = idx
	}
}

func columnMeta(h Headers, n int) (ColumnMeta, error) {
	key := func(prefix string) string {
		return fmt.Sprintf("%s%d", prefix, n)
	}

	format := h[key("TFORM")]
	if format == "" {
		return ColumnMeta{}, fmt.Errorf("missing TFORM%d", n)
	}

	dim, err := parseDim(h[key("TDIM")])
	if err != nil {
		return ColumnMeta{}, fmt.Errorf("invalid TDIM%d: %w", n, err)
	}

	return ColumnMeta{
		Index:   n - 1,
		Name:    h[key("TTYPE")],
		Format:  format,
		Unit:    h[key("TUNIT")],
		Dim:     dim,
		Scale:   floatHeader(h, key("TSCAL"), 1),
		Zero:    floatHeader(h, key("TZERO"), 0),
		NullVal: h[key("TNULL")],
	}, nil
}

func parseDim(s string) ([]int, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, nil
	}

	if !strings.HasPrefix(s, "(") || !strings.HasSuffix(s, ")") {
		return nil, fmt.Errorf("%q", s)
	}

	parts := strings.Split(s[1:len(s)-1], ",")
	dim := make([]int, len(parts))
	for i, p := range parts {
		v, err := strconv.Atoi(strings.TrimSpace(p))
		if err != nil || v <= 0 {
			return nil, fmt.Errorf("%q", s)
		}

		dim[i] = v
	}

	return dim, nil
}

func parseAsciiFormat(s string) (code byte, width int, err error) {
	s = strings.TrimSpace(s)
	if len(s) < 2 {
		return 0, 0, fmt.Errorf("invalid TFORM %q", s)
	}

	code = s[0]
	switch code {
	case 'A', 'I', 'F', 'E', 'D':
	default:
		return 0, 0, fmt.Errorf("unsupported ASCII table format %q", s)
	}

	widthStr := s[1:]
	if dot := strings.IndexByte(widthStr, '.'); dot >= 0 {
		if code == 'A' || code == 'I' {
			return 0, 0, fmt.Errorf("invalid TFORM %q", s)
		}

		for _, ch := range widthStr[dot+1:] {
			if ch < '0' || ch > '9' {
				return 0, 0, fmt.Errorf("invalid TFORM %q", s)
			}
		}

		widthStr = widthStr[:dot]
	}

	width, err = strconv.Atoi(widthStr)
	if err != nil || width <= 0 {
		return 0, 0, fmt.Errorf("invalid TFORM %q", s)
	}

	return code, width, nil
}

func asciiColumn(data []byte, rowWidth, rowCount, off, width int, code byte, meta ColumnMeta) (Column, error) {
	field := func(row int) string {
		start := row*rowWidth + off
		return string(data[start : start+width])
	}

	switch code {
	case 'A':
		c := &TypedColumn[string]{meta: meta, Data: make([]string, rowCount)}
		for r := range rowCount {
			c.Data[r] = strings.TrimRight(field(r), " ")
		}

		return c, nil
	case 'I':
		nullVal := strings.TrimSpace(meta.NullVal)
		c := &TypedColumn[int64]{meta: meta, Data: make([]int64, rowCount), Null: make([]bool, rowCount)}

		for r := range rowCount {
			s := strings.TrimSpace(field(r))
			if s == "" || (nullVal != "" && s == nullVal) {
				c.Null[r] = true
				continue
			}

			v, err := strconv.ParseInt(s, 10, 64)
			if err != nil {
				return nil, fmt.Errorf("row %d: invalid integer %q", r, s)
			}

			c.Data[r] = v
		}
		return c, nil
	default: // 'F', 'E', 'D'
		c := &TypedColumn[float64]{meta: meta, Data: make([]float64, rowCount), Null: make([]bool, rowCount)}

		for r := range rowCount {
			s := strings.TrimSpace(field(r))
			if s == "" {
				c.Null[r] = true
				continue
			}

			// Accept the Fortran 'D' exponent separator (e.g. 1.5D+03).
			s = strings.ReplaceAll(
				strings.ReplaceAll(s, "D", "E"), "d", "e",
			)

			v, err := strconv.ParseFloat(s, 64)
			if err != nil {
				return nil, fmt.Errorf("row %d: invalid float %q", r, s)
			}

			c.Data[r] = v
		}

		return c, nil
	}
}

func parseBinFormat(s string) (repeat int, code, elem byte, err error) {
	s = strings.TrimSpace(s)

	i := 0
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
	}

	repeat = 1
	if i > 0 {
		repeat, err = strconv.Atoi(s[:i])
		if err != nil {
			return 0, 0, 0, fmt.Errorf("invalid TFORM %q", s)
		}
	}

	if repeat < 1 || i == len(s) {
		return 0, 0, 0, fmt.Errorf("invalid TFORM %q", s)
	}

	code = s[i]
	i++
	switch code {
	case 'L', 'B', 'I', 'J', 'K', 'A', 'E', 'D', 'C', 'M', 'X':
		if i != len(s) {
			return 0, 0, 0, fmt.Errorf("invalid TFORM %q", s)
		}
	case 'P', 'Q':
		if i == len(s) {
			return 0, 0, 0, fmt.Errorf("invalid TFORM %q: missing array element type", s)
		}

		elem = s[i]
		i++
		if i < len(s) {
			// Optional "(maxelem)" suffix; the maximum is not used.
			rest := s[i:]
			if !strings.HasPrefix(rest, "(") || !strings.HasSuffix(rest, ")") {
				return 0, 0, 0, fmt.Errorf("invalid TFORM %q", s)
			}
		}
	default:
		return 0, 0, 0, fmt.Errorf("unsupported BINTABLE format %q", s)
	}

	return repeat, code, elem, nil
}

func binElemSize(code byte) int {
	switch code {
	case 'L', 'B', 'A':
		return 1
	case 'I':
		return 2
	case 'J', 'E':
		return 4
	case 'K', 'D', 'C':
		return 8
	case 'M':
		return 16
	}
	return 0
}

func binColumnWidth(repeat int, code byte) (int, error) {
	switch code {
	case 'X':
		return (repeat + 7) / 8, nil
	case 'P':
		return 8, nil
	case 'Q':
		return 16, nil
	}

	es := binElemSize(code)
	if es == 0 {
		return 0, fmt.Errorf("unsupported BINTABLE format code %q", code)
	}

	return repeat * es, nil
}

func intNull(meta ColumnMeta) (func(int64) bool, error) {
	if meta.NullVal == "" {
		return func(int64) bool { return false }, nil
	}

	v, err := strconv.ParseInt(meta.NullVal, 10, 64)
	if err != nil {
		return nil, fmt.Errorf("invalid TNULL %q", meta.NullVal)
	}

	return func(x int64) bool { return x == v }, nil
}

func readInt16(b []byte) int16 {
	return int16(binary.BigEndian.Uint16(b))
}

func readInt32(b []byte) int32 {
	return int32(binary.BigEndian.Uint32(b))
}

func readInt64(b []byte) int64 {
	return int64(binary.BigEndian.Uint64(b))
}

func readFloat32(b []byte) float32 {
	return math.Float32frombits(binary.BigEndian.Uint32(b))
}

func readFloat64(b []byte) float64 {
	return math.Float64frombits(binary.BigEndian.Uint64(b))
}

func readComplex64(b []byte) complex64 {
	return complex(readFloat32(b), readFloat32(b[4:]))
}

func readComplex128(b []byte) complex128 {
	return complex(readFloat64(b), readFloat64(b[8:]))
}

func binColumn(
	u HDU,
	meta ColumnMeta,
	rowWidth, rowCount,
	off, repeat int,
	code,
	elem byte,
	heapStart int,
) (Column, error) {
	data := u.Data

	switch code {
	case 'L':
		c := &TypedColumn[bool]{meta: meta, Data: make([]bool, rowCount), Null: make([]bool, rowCount)}

		for r := range rowCount {
			switch b := data[r*rowWidth+off]; b {
			case 'T':
				c.Data[r] = true
			case 'F':
				c.Data[r] = false
			default:
				c.Null[r] = true
			}
		}

		return c, nil
	case 'A':
		c := &TypedColumn[string]{meta: meta, Data: make([]string, rowCount)}

		for r := range rowCount {
			start := r*rowWidth + off
			c.Data[r] = strings.TrimRight(string(data[start:start+repeat]), " ")
		}

		return c, nil
	case 'X':
		nbytes := (repeat + 7) / 8
		c := &TypedColumn[BitArray]{meta: meta, Data: make([]BitArray, rowCount)}

		for r := range rowCount {
			start := r*rowWidth + off
			c.Data[r] = BitArray{BitCount: repeat, Bytes: bytes.Clone(data[start : start+nbytes])}
		}

		return c, nil
	case 'B':
		null, err := intNull(meta)
		if err != nil {
			return nil, err
		}

		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 1,
			func(b []byte) uint8 { return b[0] },
			func(v uint8) bool { return null(int64(v)) }), nil
	case 'I':
		null, err := intNull(meta)
		if err != nil {
			return nil, err
		}

		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 2,
			readInt16, func(v int16) bool { return null(int64(v)) }), nil
	case 'J':
		null, err := intNull(meta)
		if err != nil {
			return nil, err
		}

		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 4,
			readInt32, func(v int32) bool { return null(int64(v)) }), nil
	case 'K':
		null, err := intNull(meta)
		if err != nil {
			return nil, err
		}

		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 8,
			readInt64, func(v int64) bool { return null(v) }), nil
	case 'E':
		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 4,
			readFloat32, func(v float32) bool { return math.IsNaN(float64(v)) }), nil
	case 'D':
		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 8,
			readFloat64, math.IsNaN), nil
	case 'C':
		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 8, readComplex64,
			func(v complex64) bool {
				return math.IsNaN(float64(real(v))) || math.IsNaN(float64(imag(v)))
			}), nil
	case 'M':
		return fixedColumn(meta, data, rowWidth, rowCount, off, repeat, 16, readComplex128,
			func(v complex128) bool { return math.IsNaN(real(v)) || math.IsNaN(imag(v)) }), nil
	case 'P', 'Q':
		if repeat != 1 {
			return nil, fmt.Errorf("repeat %d is not supported for %c columns", repeat, code)
		}

		wide := code == 'Q'
		switch elem {
		case 'L':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 1,
				func(b []byte) bool { return b[0] == 'T' })
		case 'B':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 1,
				func(b []byte) uint8 { return b[0] })
		case 'I':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 2, readInt16)
		case 'J':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 4, readInt32)
		case 'K':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 8, readInt64)
		case 'E':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 4, readFloat32)
		case 'D':
			return heapColumn(meta, data, rowWidth, rowCount, off, heapStart, wide, 8, readFloat64)
		default:
			return nil, fmt.Errorf("unsupported variable-length array element type %q", elem)
		}
	}

	return nil, fmt.Errorf("unsupported BINTABLE format code %q", code)
}

func fixedColumn[T any](
	meta ColumnMeta, data []byte, rowWidth, rowCount, off, repeat, elemSize int,
	read func([]byte) T, isNull func(T) bool,
) Column {
	if repeat == 1 {
		c := &TypedColumn[T]{meta: meta, Data: make([]T, rowCount), Null: make([]bool, rowCount)}

		for r := range rowCount {
			v := read(data[r*rowWidth+off:])
			c.Data[r] = v
			c.Null[r] = isNull(v)
		}

		return c
	}

	c := &TypedColumn[[]T]{meta: meta, Data: make([][]T, rowCount)}

	for r := range rowCount {
		row := data[r*rowWidth+off:]
		vals := make([]T, repeat)

		for j := range vals {
			vals[j] = read(row[j*elemSize:])
		}

		c.Data[r] = vals
	}

	return c
}

func heapColumn[T any](
	meta ColumnMeta, data []byte, rowWidth, rowCount, off, heapStart int,
	wide bool, elemSize int, read func([]byte) T,
) (Column, error) {
	heapLen := int64(len(data) - heapStart)
	c := &TypedColumn[[]T]{meta: meta, Data: make([][]T, rowCount)}

	for r := range rowCount {
		desc := data[r*rowWidth+off:]

		var count, offset int64
		if wide {
			count = readInt64(desc)
			offset = readInt64(desc[8:])
		} else {
			count = int64(readInt32(desc))
			offset = int64(readInt32(desc[4:]))
		}

		if count < 0 || offset < 0 || count > heapLen || offset > heapLen ||
			count*int64(elemSize) > heapLen-offset {
			return nil, fmt.Errorf(
				"row %d: array descriptor (count=%d, offset=%d) out of heap bounds",
				r, count, offset,
			)
		}

		vals := make([]T, count)
		start := heapStart + int(offset)
		for j := range vals {
			vals[j] = read(data[start+j*elemSize:])
		}

		c.Data[r] = vals
	}

	return c, nil
}
