package fits

import (
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Sizes defined by the FITS standard.
const (
	// BLOCK_SIZE is the length in bytes of a FITS logical record block.
	// Headers and data are always padded to a multiple of this size.
	BLOCK_SIZE = 2880
	// HDU_LINE_LENGTH is the length in bytes of a single header card (line).
	// Each block holds BLOCK_SIZE / HDU_LINE_LENGTH cards.
	HDU_LINE_LENGTH = 80
)

const (
	commentPrefix = "COMMENT "
	historyPrefix = "HISTORY "
	keyLength     = 8
)

const (
	maxNaxis    = 999
	maxDataSize = 1 << 60
)

var validBitpix = map[int]bool{8: true, 16: true, 32: true, -32: true, -64: true}

// ErrEmptyFile is returned when the input contains no HDUs at all,
// e.g. when the reader is empty or starts with a blank block.
var ErrEmptyFile = errors.New("failed to decode FITS file: empty file")

// Decode reads the whole FITS file from r and returns all of its HDUs in order.
// The first HDU must be a primary HDU (SIMPLE), otherwise an error is returned.
// If r contains no HDUs, ErrEmptyFile is returned. Errors from a malformed HDU
// are wrapped with its 0-based index. Data of each HDU is converted to Little Endian.
func Decode(r io.Reader) (FitsFile, error) {
	fitsFile := FitsFile{}

	for i := 0; ; i++ {
		hdu, err := decodeHDU(r)
		switch {
		case errors.Is(err, io.EOF):
			if len(fitsFile) == 0 {
				return nil, ErrEmptyFile
			}
			return fitsFile, nil
		case err != nil:
			return nil, fmt.Errorf("failed to decode HDU #%d: %w", i, err)
		default:
			if i == 0 && !hdu.IsPrimary() {
				return nil, errors.New("failed to decode FITS file: first HDU is not primary")
			}
			fitsFile = append(fitsFile, hdu)
		}
	}
}

// Decoder reads HDUs from a FITS stream one at a time. Unlike Decode, it does not
// buffer the whole file, so it is suitable for streams with many HDUs.
type Decoder struct {
	r     io.Reader
	hduIx int
}

// NewDecoder returns a Decoder that reads FITS HDUs from r.
func NewDecoder(r io.Reader) *Decoder {
	return &Decoder{r: r}
}

// HDU decodes and returns the next HDU from the stream.
// When the stream ends before any HDU was read, it returns ErrEmptyFile.
// When the stream ends after at least one HDU was read, it returns io.EOF
// (together with a zero-value HDU). Decoding errors are wrapped with the
// 0-based index of the HDU being read. The primary HDU check is not performed.
func (d *Decoder) HDU() (HDU, error) {
	hdu, err := decodeHDU(d.r)
	switch {
	case errors.Is(err, io.EOF) && d.hduIx == 0:
		return HDU{}, ErrEmptyFile
	case err != nil:
		return HDU{}, fmt.Errorf("failed to decode HDU #%d: %s", d.hduIx, err)
	}

	d.hduIx++

	return hdu, err
}

func decodeHDU(r io.Reader) (hdu HDU, err error) {
	hdu.Headers, err = decodeHeaders(r)
	if err != nil {
		return
	}

	hdu.Data, err = decodeBody(r, hdu.Headers)
	if errors.Is(err, io.EOF) {
		err = io.ErrUnexpectedEOF
	}
	return
}

func decodeHeaders(r io.Reader) (Headers, error) {
	headers := make(Headers)
	buf := make([]byte, BLOCK_SIZE)
	firstKey := ""

	for {
		if _, err := io.ReadFull(r, buf); err != nil {
			if errors.Is(err, io.EOF) && firstKey != "" {
				return nil, io.ErrUnexpectedEOF
			}
			return nil, err
		}

		if firstKey == "" && isBlank(buf) {
			return nil, io.EOF
		}

		end := false
		for i := range BLOCK_SIZE / HDU_LINE_LENGTH {
			line := string(buf[i*HDU_LINE_LENGTH : i*HDU_LINE_LENGTH+HDU_LINE_LENGTH])
			key := strings.TrimSpace(line[:keyLength])

			switch {
			case key == "END":
				end = true
			case key == "", strings.HasPrefix(line, commentPrefix), strings.HasPrefix(line, historyPrefix):
				continue
			case firstKey == "":
				firstKey = key
			}

			if end || key == "" || line[8] != '=' {
				continue
			}

			headers[key] = decodeValue(splitValue(line))
		}

		if end {
			break
		}
	}

	if firstKey != "SIMPLE" && firstKey != "XTENSION" {
		return nil, fmt.Errorf("invalid HDU: first card must be SIMPLE or XTENSION, got %q", firstKey)
	}

	return headers, nil
}

func splitValue(line string) string {
	rest := line[10:]
	inStr := false

	for i := 0; i < len(rest); i++ {
		switch {
		case rest[i] == '\'':
			inStr = !inStr
		case !inStr && rest[i] == '/':
			return rest[:i]
		}
	}

	return rest
}

func decodeValue(raw string) string {
	val := strings.TrimSpace(raw)
	if !strings.HasPrefix(val, "'") {
		return val
	}

	var sb strings.Builder
	for i := 1; i < len(val); i++ {
		if val[i] != '\'' {
			sb.WriteByte(val[i])
			continue
		}

		if i+1 < len(val) && val[i+1] == '\'' {
			sb.WriteByte('\'')
			i++
			continue
		}

		break
	}

	return strings.TrimRight(sb.String(), " ")
}

func decodeBody(r io.Reader, h Headers) ([]byte, error) {
	size, err := dataSize(h)
	if err != nil {
		return nil, err
	}

	data := make([]byte, size)
	buf := make([]byte, BLOCK_SIZE)
	off := 0

	for off < size {
		if _, err := io.ReadFull(r, buf); err != nil {
			return nil, err
		}

		n := min(BLOCK_SIZE, size-off)
		copy(data[off:], buf[:n])
		off += n
	}

	bp := bitpix(h)
	if bp < 0 {
		bp = -bp
	}
	if bs := bp / 8; bs > 1 {
		swapEndian(data, bs)
	}

	return data, nil
}

func swapEndian(data []byte, size int) {
	for i := 0; i+size <= len(data); i += size {
		for a, b := i, i+size-1; a < b; a, b = a+1, b-1 {
			data[a], data[b] = data[b], data[a]
		}
	}
}

func isBlank(b []byte) bool {
	for _, c := range b {
		if c != 0 && c != ' ' {
			return false
		}
	}
	return true
}

func intHeader(h Headers, key string) (int, error) {
	raw, ok := h[key]
	if !ok {
		return 0, fmt.Errorf("missing %s", key)
	}
	v, err := strconv.Atoi(raw)
	if err != nil {
		return 0, fmt.Errorf("invalid %s: %q", key, raw)
	}
	return v, nil
}

func dataSize(h Headers) (int, error) {
	naxis, err := intHeader(h, "NAXIS")
	if err != nil {
		return 0, err
	}

	if naxis < 0 || naxis > maxNaxis {
		return 0, fmt.Errorf("invalid NAXIS: %d", naxis)
	}

	if naxis == 0 {
		return 0, nil
	}

	bp, err := intHeader(h, "BITPIX")
	if err != nil {
		return 0, err
	}

	if !validBitpix[bp] {
		return 0, fmt.Errorf("invalid BITPIX: %d", bp)
	}

	size := int64(1)
	for i := 1; i <= naxis; i++ {
		dim, err := intHeader(h, fmt.Sprintf("NAXIS%d", i))
		if err != nil {
			return 0, err
		}
		if dim <= 0 {
			return 0, fmt.Errorf("invalid NAXIS%d: %d", i, dim)
		}
		if size > maxDataSize/int64(dim) {
			return 0, fmt.Errorf("data size exceeds maximum %d", maxDataSize)
		}
		size *= int64(dim)
	}

	if bp < 0 {
		bp = -bp
	}

	size *= int64(bp / 8)
	if size > maxDataSize {
		return 0, fmt.Errorf("data size %d exceeds maximum %d", size, maxDataSize)
	}

	return int(size), nil
}
