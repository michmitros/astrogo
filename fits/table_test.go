package fits

import (
	"encoding/binary"
	"math"
	"slices"
	"strconv"
	"strings"
	"testing"
)

func TestToTable_Acii(t *testing.T) {
	tests := map[string]struct {
		hdu     HDU
		wantErr bool
		check   func(t *testing.T, tbl *Table)
	}{
		"ascii string": {
			hdu: asciiTableHDU(5, 3, []asciiColSpec{
				{form: "A5", bcol: 1, values: []string{"abc  ", "de   ", "fgh  "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				want := []string{"abc", "de", "fgh"}
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[string]).Data, want) {
					t.Errorf("string=%v", tbl.Columns[0].(*TypedColumn[string]).Data)
				}
			},
		},
		"ascii int": {
			hdu: asciiTableHDU(5, 3, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 ", "  -7 ", "  0  "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[int64]).Data, []int64{42, -7, 0}) {
					t.Error("int mismatch")
				}
			},
		},
		"ascii float F": {
			hdu: asciiTableHDU(6, 2, []asciiColSpec{
				{form: "F6.1", bcol: 1, values: []string{"  3.14", " -1.5 "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float64])
				if math.Abs(c.Data[0]-3.14) > 0.001 || math.Abs(c.Data[1]+1.5) > 0.001 {
					t.Errorf("float=%v", c.Data)
				}
			},
		},
		"ascii float D exp": {
			hdu: asciiTableHDU(8, 2, []asciiColSpec{
				{form: "E8.2", bcol: 1, values: []string{" 1.5D+03", " 2.0d-01"}},
			}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float64])
				if math.Abs(c.Data[0]-1500) > 0.1 || math.Abs(c.Data[1]-0.2) > 0.001 {
					t.Errorf("D float=%v", c.Data)
				}
			},
		},
		"ascii int null blank": {
			hdu: asciiTableHDU(5, 3, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 ", "     ", "  0  "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[int64])
				if !c.Null[1] || c.Null[0] || c.Null[2] {
					t.Errorf("Null=%v want [false true false]", c.Null)
				}
			},
		},
		"ascii int null TNULL": {
			hdu: asciiTableHDUExtra(5, 3, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  99 ", "  42 ", "  99 "}},
			}, map[int]map[string]string{1: {"TNULL": "99"}}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[int64])
				if !c.Null[0] || !c.Null[2] || c.Null[1] || c.Data[1] != 42 {
					t.Errorf("TNULL: Null=%v Data[1]=%d", c.Null, c.Data[1])
				}
			},
		},
		"ascii float null blank": {
			hdu: asciiTableHDU(6, 3, []asciiColSpec{
				{form: "F6.1", bcol: 1, values: []string{"  3.14", "      ", "  2.0 "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				if !tbl.Columns[0].(*TypedColumn[float64]).Null[1] {
					t.Error("row 1 should be null")
				}
			},
		},
		"ascii multi col": {
			hdu: asciiTableHDU(11, 2, []asciiColSpec{
				{form: "A3", bcol: 1, values: []string{"abc", "xyz"}},
				{form: "I4", bcol: 4, values: []string{"  42", "  -1"}},
				{form: "F4.1", bcol: 8, values: []string{" 1.5", " 2.0"}},
			}),
			check: func(t *testing.T, tbl *Table) {
				if len(tbl.Columns) != 3 {
					t.Fatal("want 3 cols")
				}
				s := tbl.Columns[0].(*TypedColumn[string]).Data
				i := tbl.Columns[1].(*TypedColumn[int64]).Data
				if !slices.Equal(s, []string{"abc", "xyz"}) || !slices.Equal(i, []int64{42, -1}) {
					t.Error("multi-col mismatch")
				}
			},
		},
		"ascii TDIM": {
			hdu: asciiTableHDUExtra(5, 1, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 "}},
			}, map[int]map[string]string{1: {"TDIM": "(2,3)"}}),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].Meta().Dim, []int{2, 3}) {
					t.Errorf("Dim=%v", tbl.Columns[0].Meta().Dim)
				}
			},
		},
		"ascii TSCAL TZERO": {
			hdu: asciiTableHDUExtra(5, 1, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 "}},
			}, map[int]map[string]string{1: {"TSCAL": "2.5", "TZERO": "10.0"}}),
			check: func(t *testing.T, tbl *Table) {
				m := tbl.Columns[0].Meta()
				if m.Scale != 2.5 || m.Zero != 10.0 {
					t.Errorf("Scale=%v Zero=%v", m.Scale, m.Zero)
				}
				if tbl.Columns[0].(*TypedColumn[int64]).Data[0] != 42 {
					t.Error("data should not be scaled")
				}
			},
		},
		"ascii EXTNAME": {
			hdu: newAsciiTableHDU("FOO", 5, 1, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 "}},
			}),
			check: func(t *testing.T, tbl *Table) {
				if tbl.ExtName != "FOO" {
					t.Errorf("ExtName=%q", tbl.ExtName)
				}
			},
		},
		"ascii zero rows": {
			hdu: asciiTableHDU(5, 0, []asciiColSpec{{form: "A5", bcol: 1}}),
			check: func(t *testing.T, tbl *Table) {
				if tbl.RowCount != 0 || tbl.Columns[0].Len() != 0 {
					t.Error("zero rows mismatch")
				}
			},
		},
		"ascii non-zero PCOUNT": {
			hdu:     asciiTableHDUWithPcount(5, 1, 10),
			wantErr: true,
		},
		"ascii missing TFORM": {
			hdu: newHDU(Headers{
				"XTENSION": "TABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "5", "NAXIS2": "1", "TFIELDS": "1", "PCOUNT": "0",
			}, dryPadData([]byte("     "))),
			wantErr: true,
		},
		"ascii unsupported format": {
			hdu:     asciiTableHDU(5, 1, []asciiColSpec{{form: "X5", bcol: 1, values: []string{"     "}}}),
			wantErr: true,
		},
		"ascii TBCOL exceeds row": {
			hdu: newHDU(Headers{
				"XTENSION": "TABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "5", "NAXIS2": "1", "TFIELDS": "1", "PCOUNT": "0",
				"TFORM1": "A3", "TBCOL1": "4",
			}, dryPadData([]byte("abcde"))),
			wantErr: true,
		},
		"ascii invalid integer": {
			hdu:     asciiTableHDU(5, 1, []asciiColSpec{{form: "I5", bcol: 1, values: []string{"abc  "}}}),
			wantErr: true,
		},
		"ascii invalid float": {
			hdu:     asciiTableHDU(6, 1, []asciiColSpec{{form: "F6.1", bcol: 1, values: []string{"abc   "}}}),
			wantErr: true,
		},
		"invalid TDIM": {
			hdu: asciiTableHDUExtra(5, 1, []asciiColSpec{
				{form: "I5", bcol: 1, values: []string{"  42 "}},
			}, map[int]map[string]string{1: {"TDIM": "2,3"}}),
			wantErr: true,
		},
		"TBCOL zero": {
			hdu: newHDU(Headers{
				"XTENSION": "TABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "5", "NAXIS2": "1", "TFIELDS": "1", "PCOUNT": "0",
				"TFORM1": "A5", "TBCOL1": "0",
			}, dryPadData([]byte("abc  "))),
			wantErr: true,
		},
		// --- Error cases ---
		"non-table HDU": {
			hdu: newHDU(Headers{
				"XTENSION": "IMAGE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "2", "NAXIS2": "2",
			}, []byte{0, 0, 0, 0}),
			wantErr: true,
		},
		"NAXIS not 2": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "8", "NAXIS": "1", "NAXIS1": "4",
			}, nil),
			wantErr: true,
		},
		"BITPIX not 8": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "16", "NAXIS": "2",
				"NAXIS1": "4", "NAXIS2": "1", "TFIELDS": "1", "PCOUNT": "0",
			}, nil),
			wantErr: true,
		},
		"negative NAXIS1": {
			hdu: binHDUExtra(4, 1, "1J", bigInt32(42),
				map[string]string{"NAXIS1": "-1"}),
			wantErr: true,
		},
		"missing TFIELDS": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "4", "NAXIS2": "1", "PCOUNT": "0",
			}, nil),
			wantErr: true,
		},
		"negative PCOUNT": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "4", "NAXIS2": "1", "TFIELDS": "1", "PCOUNT": "-1",
			}, nil),
			wantErr: true,
		},
		"data too short heap": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "4", "NAXIS2": "1", "TFIELDS": "1",
				"PCOUNT": "100", "TFORM1": "1J",
			}, bigInt32(42)),
			wantErr: true,
		},
		"col width exceeds row": {
			hdu: newHDU(Headers{
				"XTENSION": "BINTABLE", "BITPIX": "8", "NAXIS": "2",
				"NAXIS1": "4", "NAXIS2": "1", "TFIELDS": "2",
				"PCOUNT": "0", "TFORM1": "1J", "TFORM2": "1J",
			}, bigInt32(1, 2)),
			wantErr: true,
		},
		"data too short": {
			hdu:     binHDU(4, 2, "1J", bigInt32(42), nil),
			wantErr: true,
		},
		"invalid TNULL": {
			hdu: binHDUExtra(1, 1, "1B", []byte{42},
				map[string]string{"TNULL": "abc"}),
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			tbl, err := ToTable(tt.hdu)
			if (err != nil) != tt.wantErr {
				t.Fatalf("ToTable() error=%v wantErr=%v", err, tt.wantErr)
			}
			if tt.wantErr {
				return
			}
			if tt.check != nil {
				tt.check(t, tbl)
			}
		})
	}
}

func TestToTable_BinTable(t *testing.T) {
	tests := map[string]struct {
		hdu     HDU
		wantErr bool
		check   func(t *testing.T, tbl *Table)
	}{
		"bin boolean L": {
			hdu: binHDU(1, 3, "1L", []byte{'T', 'F', 'X'}, nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[bool])
				if !c.Data[0] || c.Data[1] || !c.Null[2] {
					t.Errorf("bool=%v null=%v", c.Data, c.Null)
				}
			},
		},
		"bin uint8 B": {
			hdu: binHDU(1, 3, "1B", []byte{0, 128, 255}, nil),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[uint8]).Data, []uint8{0, 128, 255}) {
					t.Error("uint8 mismatch")
				}
			},
		},
		"bin int16 I": {
			hdu: binHDU(2, 3, "1I", bigInt16(0, 256, -32768), nil),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[int16]).Data, []int16{0, 256, -32768}) {
					t.Error("int16 mismatch")
				}
			},
		},
		"bin int32 J": {
			hdu: binHDU(4, 2, "1J", bigInt32(42, -1), nil),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[int32]).Data, []int32{42, -1}) {
					t.Error("int32 mismatch")
				}
			},
		},
		"bin int64 K": {
			hdu: binHDU(8, 2, "1K", bigInt64(1, -9999), nil),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[int64]).Data, []int64{1, -9999}) {
					t.Error("int64 mismatch")
				}
			},
		},
		"bin float32 E": {
			hdu: binHDU(4, 2, "1E", bigFloat32(3.14, -2.5), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float32])
				if math.Abs(float64(c.Data[0])-3.14) > 0.001 || math.Abs(float64(c.Data[1])+2.5) > 0.001 {
					t.Errorf("float32=%v", c.Data)
				}
			},
		},
		"bin float64 D": {
			hdu: binHDU(8, 2, "1D", bigFloat64(1.5, -0.25), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float64])
				if math.Abs(c.Data[0]-1.5) > 1e-9 || math.Abs(c.Data[1]+0.25) > 1e-9 {
					t.Errorf("float64=%v", c.Data)
				}
			},
		},
		"bin complex64 C": {
			hdu: binHDU(8, 1, "1C", bigComplex64(complex64(complex(1.5, 2.5))), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[complex64])
				if c.Data[0] != complex64(complex(1.5, 2.5)) {
					t.Errorf("complex64=%v", c.Data[0])
				}
			},
		},
		"bin complex128 M": {
			hdu: binHDU(16, 1, "1M", bigComplex128(complex(1.5, -3.25)), nil),
			check: func(t *testing.T, tbl *Table) {
				if tbl.Columns[0].(*TypedColumn[complex128]).Data[0] != complex(1.5, -3.25) {
					t.Error("complex128 mismatch")
				}
			},
		},
		"bin string A": {
			hdu: binHDU(5, 2, "5A", []byte("ab cdxy z  "), nil),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[string]).Data, []string{"ab cd", "xy z"}) {
					t.Error("string mismatch")
				}
			},
		},
		"bin bit array X": {
			hdu: binHDU(1, 2, "5X", []byte{0b10100000, 0b01000000}, nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[BitArray])
				if c.Data[0].BitCount != 5 || c.Data[0].Bytes[0] != 0b10100000 {
					t.Errorf("BitArray=%+v", c.Data[0])
				}
			},
		},
		"bin multi int16 r=3": {
			hdu: binHDU(6, 2, "3I", bigInt16(1, 2, 3, 4, 5, 6), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[[]int16])
				w := [][]int16{{1, 2, 3}, {4, 5, 6}}
				for i := range w {
					if !slices.Equal(c.Data[i], w[i]) {
						t.Errorf("row %d=%v", i, c.Data[i])
					}
				}
			},
		},
		"bin multi uint8 r=2": {
			hdu: binHDU(2, 2, "2B", []byte{10, 20, 30, 40}, nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[[]uint8])
				w := [][]uint8{{10, 20}, {30, 40}}
				for i := range w {
					if !slices.Equal(c.Data[i], w[i]) {
						t.Errorf("row %d=%v", i, c.Data[i])
					}
				}
			},
		},
		// nulls
		"bin null float32 NaN": {
			hdu: binHDU(4, 3, "1E", bigFloat32(1.0, float32(math.NaN()), -2.0), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float32])
				if !c.Null[1] || c.Null[0] || c.Null[2] {
					t.Errorf("Null=%v want [false true false]", c.Null)
				}
			},
		},
		"bin null float64 NaN": {
			hdu: binHDU(8, 2, "1D", bigFloat64(math.NaN(), 42.0), nil),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[float64])
				if !c.Null[0] || c.Null[1] {
					t.Errorf("Null=%v want [true false]", c.Null)
				}
			},
		},
		"bin null complex64 NaN": {
			hdu: binHDU(8, 1, "1C", bigComplex64(complex(float32(math.NaN()), 1.0)), nil),
			check: func(t *testing.T, tbl *Table) {
				if !tbl.Columns[0].(*TypedColumn[complex64]).Null[0] {
					t.Error("should be null")
				}
			},
		},
		"bin null complex128 NaN": {
			hdu: binHDU(16, 1, "1M", bigComplex128(complex(1.0, math.NaN())), nil),
			check: func(t *testing.T, tbl *Table) {
				if !tbl.Columns[0].(*TypedColumn[complex128]).Null[0] {
					t.Error("should be null")
				}
			},
		},
		"bin TNULL uint8": {
			hdu: binHDUExtra(1, 3, "1B", []byte{99, 42, 99}, map[string]string{"TNULL": "99"}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[uint8])
				if !c.Null[0] || c.Null[1] || !c.Null[2] || c.Data[1] != 42 {
					t.Errorf("Null=%v Data[1]=%d", c.Null, c.Data[1])
				}
			},
		},
		"bin TNULL int16": {
			hdu: binHDUExtra(2, 2, "1I", bigInt16(-99, 42), map[string]string{"TNULL": "-99"}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[int16])
				if !c.Null[0] || c.Null[1] {
					t.Errorf("Null=%v", c.Null)
				}
			},
		},
		"bin TNULL int32": {
			hdu: binHDUExtra(4, 2, "1J", bigInt32(0, -1), map[string]string{"TNULL": "-1"}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[int32])
				if c.Null[0] || !c.Null[1] {
					t.Errorf("Null=%v", c.Null)
				}
			},
		},
		"bin TNULL int64": {
			hdu: binHDUExtra(8, 2, "1K", bigInt64(999, 999), map[string]string{"TNULL": "999"}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[int64])
				if !c.Null[0] || !c.Null[1] {
					t.Errorf("Null=%v", c.Null)
				}
			},
		},
		"bin empty TNULL": {
			hdu: binHDUExtra(1, 2, "1B", []byte{42, 99}, map[string]string{"TNULL": ""}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[uint8])
				if c.Null[0] || c.Null[1] {
					t.Errorf("no nulls: %v", c.Null)
				}
			},
		},
		// var length
		"bin var P uint8": {
			hdu: binVarHDU(1, false, 'B', [][]byte{{10, 20, 30}}),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[[]uint8]).Data[0], []uint8{10, 20, 30}) {
					t.Error("heap uint8 mismatch")
				}
			},
		},
		"bin var P int16": {
			hdu: binVarHDU(2, false, 'I', [][]byte{bigInt16(1, -2)}),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[[]int16]).Data[0], []int16{1, -2}) {
					t.Error("heap int16 mismatch")
				}
			},
		},
		"bin var Q float64": {
			hdu: binVarHDU(8, true, 'D', [][]byte{bigFloat64(1.5), bigFloat64(-2.5, 3.0)}),
			check: func(t *testing.T, tbl *Table) {
				c := tbl.Columns[0].(*TypedColumn[[]float64])
				if len(c.Data[0]) != 1 || math.Abs(c.Data[0][0]-1.5) > 1e-9 {
					t.Errorf("row0=%v", c.Data[0])
				}
				if len(c.Data[1]) != 2 || math.Abs(c.Data[1][1]-3.0) > 1e-9 {
					t.Errorf("row1=%v", c.Data[1])
				}
			},
		},
		"bin var P bool": {
			hdu: binVarHDU(1, false, 'L', [][]byte{{'T', 'F', 'T'}}),
			check: func(t *testing.T, tbl *Table) {
				if !slices.Equal(tbl.Columns[0].(*TypedColumn[[]bool]).Data[0], []bool{true, false, true}) {
					t.Error("heap bool mismatch")
				}
			},
		},
		"bin var repeat>1 error": {
			hdu:     binHDUExtra(8, 1, "2PB", make([]byte, 8), nil),
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			tbl, err := ToTable(tt.hdu)
			if (err != nil) != tt.wantErr {
				t.Fatalf("ToTable() error=%v wantErr=%v", err, tt.wantErr)
			}
			if tt.wantErr {
				return
			}
			if tt.check != nil {
				tt.check(t, tbl)
			}
		})
	}
}

func TestGetColumn(t *testing.T) {
	rows := slices.Concat(bigInt32(1), bigFloat64(3.0), bigInt32(2), bigFloat64(4.0))
	tbl, _ := ToTable(binHDUExtra(12, 2, "1J,1D", rows,
		map[string]string{"TTYPE1": "ALPHA", "TTYPE2": "BETA"}))

	t.Run("get int32", func(t *testing.T) {
		data, err := GetColumn[int32](tbl, "ALPHA")
		if err != nil || !slices.Equal(data, []int32{1, 2}) {
			t.Errorf("GetColumn ALPHA: err=%v data=%v", err, data)
		}
	})
	t.Run("get float64", func(t *testing.T) {
		data, err := GetColumn[float64](tbl, "BETA")
		if err != nil || !slices.Equal(data, []float64{3.0, 4.0}) {
			t.Errorf("GetColumn BETA: err=%v data=%v", err, data)
		}
	})
	t.Run("missing", func(t *testing.T) {
		if _, err := GetColumn[int32](tbl, "NOPE"); err == nil {
			t.Error("should error")
		}
	})
	t.Run("wrong type", func(t *testing.T) {
		if _, err := GetColumn[float64](tbl, "ALPHA"); err == nil {
			t.Error("should error")
		}
	})
}

func TestTypedColumn(t *testing.T) {
	meta := ColumnMeta{Index: 2, Name: "FOO", Format: "1J", Unit: "m", Dim: []int{2, 3}, Scale: 2.0, Zero: 10.0, NullVal: "-99"}
	col := &TypedColumn[int32]{meta: meta, Data: []int32{1, 2, 3}}

	t.Run("Meta", func(t *testing.T) {
		m := col.Meta()
		if m.Index != 2 || m.Name != "FOO" || m.Format != "1J" || m.Unit != "m" || !slices.Equal(m.Dim, []int{2, 3}) || m.Scale != 2.0 || m.Zero != 10.0 || m.NullVal != "-99" {
			t.Errorf("Meta=%+v", m)
		}
	})
	t.Run("Len", func(t *testing.T) {
		if col.Len() != 3 {
			t.Errorf("Len=%d", col.Len())
		}
	})
	t.Run("GetAny value", func(t *testing.T) {
		if col.GetAny(0) != int32(1) {
			t.Errorf("GetAny(0)=%v", col.GetAny(0))
		}
	})
	t.Run("GetAny nil", func(t *testing.T) {
		nc := &TypedColumn[float64]{Data: []float64{1.0, 2.0}, Null: []bool{true, false}}
		if nc.GetAny(0) != nil {
			t.Errorf("GetAny(0)=%v want nil", nc.GetAny(0))
		}
		if nc.GetAny(1) != 2.0 {
			t.Errorf("GetAny(1)=%v", nc.GetAny(1))
		}
	})
	t.Run("GetAny no Null", func(t *testing.T) {
		nc := &TypedColumn[string]{Data: []string{"hello"}}
		if nc.GetAny(0) != "hello" {
			t.Errorf("GetAny(0)=%v", nc.GetAny(0))
		}
	})
}

func TestParseDim(t *testing.T) {
	tests := map[string]struct {
		s       string
		want    []int
		wantErr bool
	}{
		"empty":         {s: "", want: nil},
		"(2,3)":         {s: "(2,3)", want: []int{2, 3}},
		"(10)":          {s: "(10)", want: []int{10}},
		"no parens":     {s: "2,3", wantErr: true},
		"non numeric":   {s: "(a,3)", wantErr: true},
		"zero":          {s: "(0,3)", wantErr: true},
		"negative":      {s: "(-1,3)", wantErr: true},
		"spaces around": {s: "  (2 , 3)  ", want: []int{2, 3}},
	}
	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := parseDim(tt.s)

			if (err != nil) != tt.wantErr {
				t.Fatalf("parseDim err=%v wantErr=%v", err, tt.wantErr)
			}

			if !slices.Equal(got, tt.want) {
				t.Errorf("parseDim=%v want=%v", got, tt.want)
			}
		})
	}
}

func TestParseAsciiFormat(t *testing.T) {
	tests := map[string]struct {
		s       string
		want    byte
		width   int
		wantErr bool
	}{
		"A10":         {s: "A10", want: 'A', width: 10},
		"I5":          {s: "I5", want: 'I', width: 5},
		"F10.2":       {s: "F10.2", want: 'F', width: 10},
		"E12.3":       {s: "E12.3", want: 'E', width: 12},
		"D20.6":       {s: "D20.6", want: 'D', width: 20},
		"too short":   {s: "A", wantErr: true},
		"unsupported": {s: "X5", wantErr: true},
		"zero width":  {s: "I0", wantErr: true},
		"negative":    {s: "I-5", wantErr: true},
		"A with dot":  {s: "A5.2", wantErr: true},
		"I with dot":  {s: "I5.2", wantErr: true},
		"non numeric": {s: "Iab", wantErr: true},
		"bad decimal": {s: "F5.xy", wantErr: true},
		"spaces":      {s: "  A10  ", want: 'A', width: 10},
	}
	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			code, width, err := parseAsciiFormat(tt.s)

			if (err != nil) != tt.wantErr {
				t.Fatalf("parseAsciiFormat err=%v wantErr=%v", err, tt.wantErr)
			}

			if tt.wantErr {
				return
			}

			if code != tt.want || width != tt.width {
				t.Errorf("parseAsciiFormat=(%c, %d) want=(%c, %d)", code, width, tt.want, tt.width)
			}
		})
	}
}

// asciiColSpec describes a single column of a test ASCII table. bcol is the
// 1-based starting column (TBCOLn) and values holds the raw field text per row.
type asciiColSpec struct {
	form   string
	bcol   int
	values []string
}

// newHDU creates an HDU from headers and raw data.
func newHDU(h Headers, data []byte) HDU {
	return HDU{Headers: h, Data: data}
}

// dryPadData returns a copy of data without any block padding (HDU.Data
// never contains padding).
func dryPadData(data []byte) []byte {
	return slices.Clone(data)
}

// newAsciiTableHDU builds an ASCII TABLE HDU with the given EXTNAME.
func newAsciiTableHDU(extName string, rowWidth, rows int, cols []asciiColSpec) HDU {
	u := asciiTableHDU(rowWidth, rows, cols)
	u.Headers["EXTNAME"] = extName
	return u
}

// asciiTableHDU builds an ASCII TABLE HDU of rows rows, each rowWidth bytes
// wide. Fields are placed at TBCOLn; the remaining bytes are spaces.
func asciiTableHDU(rowWidth, rows int, cols []asciiColSpec) HDU {
	return asciiTableHDUExtra(rowWidth, rows, cols, nil)
}

// asciiTableHDUExtra works like asciiTableHDU and additionally sets extra
// per-column keywords: extra[n][key] becomes the header key+n (e.g. TNULL1).
func asciiTableHDUExtra(rowWidth, rows int, cols []asciiColSpec, extra map[int]map[string]string) HDU {
	h := Headers{
		"XTENSION": "TABLE",
		"BITPIX":   "8",
		"NAXIS":    "2",
		"NAXIS1":   strconv.Itoa(rowWidth),
		"NAXIS2":   strconv.Itoa(rows),
		"PCOUNT":   "0",
		"GCOUNT":   "1",
		"TFIELDS":  strconv.Itoa(len(cols)),
	}

	data := []byte(strings.Repeat(" ", rowWidth*rows))
	for i, c := range cols {
		n := strconv.Itoa(i + 1)
		h["TFORM"+n] = c.form
		h["TBCOL"+n] = strconv.Itoa(c.bcol)

		for r, v := range c.values {
			if r < rows && c.bcol >= 1 {
				copy(data[r*rowWidth+c.bcol-1:], v)
			}
		}
	}

	for n, kv := range extra {
		for k, v := range kv {
			h[k+strconv.Itoa(n)] = v
		}
	}

	return newHDU(h, data)
}

// asciiTableHDUWithPcount builds a single-column ASCII TABLE HDU with the
// given (invalid) PCOUNT.
func asciiTableHDUWithPcount(rowWidth, rows, pcount int) HDU {
	u := asciiTableHDU(rowWidth, rows, []asciiColSpec{
		{form: "A" + strconv.Itoa(rowWidth), bcol: 1},
	})

	u.Headers["PCOUNT"] = strconv.Itoa(pcount)
	return u
}

// binHDU builds a BINTABLE HDU with the given row width in bytes, row count,
// comma-separated TFORMs and raw data.
func binHDU(rowWidth, rows int, tform string, data []byte, extra map[string]string) HDU {
	return binHDUExtra(rowWidth, rows, tform, data, extra)
}

// binHDUExtra works like binHDU and additionally sets extra header keywords.
// Keys ending with a digit (e.g. TTYPE2, NAXIS1) are used as is; other keys
// (e.g. TNULL) are applied to the first column (TNULL1).
func binHDUExtra(rowWidth, rows int, tform string, data []byte, extra map[string]string) HDU {
	forms := strings.Split(tform, ",")
	h := Headers{
		"XTENSION": "BINTABLE",
		"BITPIX":   "8",
		"NAXIS":    "2",
		"NAXIS1":   strconv.Itoa(rowWidth),
		"NAXIS2":   strconv.Itoa(rows),
		"PCOUNT":   "0",
		"GCOUNT":   "1",
		"TFIELDS":  strconv.Itoa(len(forms)),
	}

	for i, f := range forms {
		h["TFORM"+strconv.Itoa(i+1)] = f
	}

	for k, v := range extra {
		if last := k[len(k)-1]; last < '0' || last > '9' {
			k += "1"
		}

		h[k] = v
	}

	return newHDU(h, data)
}

// binVarHDU builds a BINTABLE HDU with a single variable-length array column
// ('P' or, if wide, 'Q') of elem type. elemSize is the element size in bytes;
// rows holds the raw heap bytes of every row's array.
func binVarHDU(elemSize int, wide bool, elem byte, rows [][]byte) HDU {
	var desc, heap []byte
	for _, r := range rows {
		count, offset := len(r)/elemSize, len(heap)
		if wide {
			desc = binary.BigEndian.AppendUint64(desc, uint64(count))
			desc = binary.BigEndian.AppendUint64(desc, uint64(offset))
		} else {
			desc = binary.BigEndian.AppendUint32(desc, uint32(count))
			desc = binary.BigEndian.AppendUint32(desc, uint32(offset))
		}

		heap = append(heap, r...)
	}

	code := "P"
	if wide {
		code = "Q"
	}

	u := binHDU(len(desc)/max(len(rows), 1), len(rows), "1"+code+string(elem), append(desc, heap...), nil)
	u.Headers["PCOUNT"] = strconv.Itoa(len(heap))

	return u
}

func bigInt16(vals ...int16) []byte {
	var b []byte
	for _, v := range vals {
		b = binary.BigEndian.AppendUint16(b, uint16(v))
	}

	return b
}

func bigInt32(vals ...int32) []byte {
	var b []byte
	for _, v := range vals {
		b = binary.BigEndian.AppendUint32(b, uint32(v))
	}

	return b
}

func bigInt64(vals ...int64) []byte {
	var b []byte
	for _, v := range vals {
		b = binary.BigEndian.AppendUint64(b, uint64(v))
	}

	return b
}

func bigFloat32(vals ...float32) []byte {
	var b []byte
	for _, v := range vals {
		b = binary.BigEndian.AppendUint32(b, math.Float32bits(v))
	}

	return b
}

func bigFloat64(vals ...float64) []byte {
	var b []byte
	for _, v := range vals {
		b = binary.BigEndian.AppendUint64(b, math.Float64bits(v))
	}

	return b
}

func bigComplex64(vals ...complex64) []byte {
	var b []byte
	for _, v := range vals {
		b = append(b, bigFloat32(real(v), imag(v))...)
	}

	return b
}

func bigComplex128(vals ...complex128) []byte {
	var b []byte
	for _, v := range vals {
		b = append(b, bigFloat64(real(v), imag(v))...)
	}

	return b
}
