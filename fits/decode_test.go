package fits

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"maps"
	"strings"
	"testing"
)

func TestDecode(t *testing.T) {
	tests := map[string]struct {
		stream    []byte
		wantData  [][]byte
		wantErr   bool
		wantErrIs error
	}{
		"empty reader": {
			stream:    nil,
			wantErr:   true,
			wantErrIs: ErrEmptyFile,
		},
		"primary HDU without data": {
			stream:   hduBytes([]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")}, nil),
			wantData: [][]byte{nil},
		},
		"primary HDU with 8-bit data": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "1"), card("NAXIS1", "3")},
				[]byte{1, 2, 3},
			),
			wantData: [][]byte{{1, 2, 3}},
		},
		"16-bit data converted to little endian": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "16"), card("NAXIS", "1"), card("NAXIS1", "2")},
				[]byte{0x01, 0x02, 0x03, 0x04},
			),
			wantData: [][]byte{{0x02, 0x01, 0x04, 0x03}},
		},
		"primary HDU and image extension": {
			stream: append(
				hduBytes([]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")}, nil),
				hduBytes([]string{
					card("XTENSION", "'IMAGE   '"),
					card("BITPIX", "8"),
					card("NAXIS", "1"),
					card("NAXIS1", "2"),
					card("PCOUNT", "0"),
					card("GCOUNT", "1"),
				}, []byte{9, 8})...,
			),
			wantData: [][]byte{nil, {9, 8}},
		},
		"first HDU is an extension": {
			stream: hduBytes(
				[]string{card("XTENSION", "'IMAGE   '"), card("BITPIX", "8"), card("NAXIS", "0")},
				nil,
			),
			wantErr: true,
		},
		"first card is neither SIMPLE nor XTENSION": {
			stream:  hduBytes([]string{raw("FOO     = T"), card("BITPIX", "8"), card("NAXIS", "0")}, nil),
			wantErr: true,
		},
		"truncated data": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "1"), card("NAXIS1", "10")},
				nil,
			),
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
		"second HDU is invalid": {
			stream: append(
				hduBytes([]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")}, nil),
				hduBytes([]string{card("XTENSION", "'IMAGE   '"), card("NAXIS", "-1")}, nil)...,
			),
			wantErr: true,
		},
		"header without END card": {
			stream:    padBlock([]byte(card("SIMPLE", "T")+card("BITPIX", "8")), ' '),
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			file, err := Decode(bytes.NewReader(tt.stream))
			if (err != nil) != tt.wantErr {
				t.Fatalf("Decode() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErrIs != nil && !errors.Is(err, tt.wantErrIs) {
				t.Fatalf("Decode() error = %v, want error matching %v", err, tt.wantErrIs)
			}
			if tt.wantErr {
				return
			}

			if len(file) != len(tt.wantData) {
				t.Fatalf("Decode() returned %d HDUs, want %d", len(file), len(tt.wantData))
			}
			for i := range file {
				if !bytes.Equal(file[i].Data, tt.wantData[i]) {
					t.Errorf("Decode() HDU #%d data = %v, want %v", i, file[i].Data, tt.wantData[i])
				}
			}
			if !file[0].IsPrimary() {
				t.Errorf("Decode() first HDU is not primary")
			}
		})
	}
}

func TestDecoderHDU(t *testing.T) {
	tests := map[string]struct {
		stream    []byte
		wantCount int
		wantErrIs error
	}{
		"empty stream": {
			stream:    nil,
			wantCount: 0,
			wantErrIs: ErrEmptyFile,
		},
		"single primary HDU": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")},
				nil,
			),
			wantCount: 1,
		},
		"two HDUs": {
			stream: append(
				hduBytes([]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")}, nil),
				hduBytes([]string{
					card("XTENSION", "'IMAGE   '"),
					card("BITPIX", "8"),
					card("NAXIS", "1"),
					card("NAXIS1", "2"),
				}, []byte{1, 2})...,
			),
			wantCount: 2,
		},
		"invalid HDU": {
			stream:    hduBytes([]string{card("SIMPLE", "T"), card("NAXIS", "-1")}, nil),
			wantCount: 0,
		},
		"non primary first HDU is accepted": {
			stream: hduBytes(
				[]string{card("XTENSION", "'IMAGE   '"), card("BITPIX", "8"), card("NAXIS", "0")},
				nil,
			),
			wantCount: 1,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			d := NewDecoder(bytes.NewReader(tt.stream))

			count := 0
			var err error
			for range 10 {
				_, err = d.HDU()
				if err != nil {
					break
				}
				count++
			}

			if err == nil {
				t.Fatalf("HDU() expected an error after %d HDUs, got none", count)
			}
			if count != tt.wantCount {
				t.Fatalf("HDU() read %d HDUs, want %d (last error: %v)", count, tt.wantCount, err)
			}
			if tt.wantErrIs != nil && !errors.Is(err, tt.wantErrIs) {
				t.Fatalf("HDU() error = %v, want error matching %v", err, tt.wantErrIs)
			}
		})
	}
}

func TestNewDecoder(t *testing.T) {
	tests := map[string]struct {
		stream []byte
	}{
		"empty stream": {
			stream: nil,
		},
		"stream with data": {
			stream: []byte("SIMPLE  = T"),
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			d := NewDecoder(bytes.NewReader(tt.stream))
			if d == nil {
				t.Fatal("NewDecoder() returned nil")
			}
			if d.hduIx != 0 {
				t.Errorf("NewDecoder() hduIx = %d, want 0", d.hduIx)
			}

			got, err := io.ReadAll(d.r)
			if err != nil {
				t.Fatalf("io.ReadAll() error = %v", err)
			}
			if !bytes.Equal(got, tt.stream) {
				t.Errorf("NewDecoder() reader content = %q, want %q", got, tt.stream)
			}
		})
	}
}

func TestDecodeHDU(t *testing.T) {
	tests := map[string]struct {
		stream    []byte
		wantData  []byte
		wantErr   bool
		wantErrIs error
	}{
		"valid 8-bit HDU": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "1"), card("NAXIS1", "2")},
				[]byte{5, 6},
			),
			wantData: []byte{5, 6},
		},
		"empty stream": {
			stream:    nil,
			wantErr:   true,
			wantErrIs: io.EOF,
		},
		"missing data": {
			stream: hduBytes(
				[]string{card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "1"), card("NAXIS1", "2")},
				nil,
			),
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
		"missing NAXIS": {
			stream:  hduBytes([]string{card("SIMPLE", "T"), card("BITPIX", "8")}, nil),
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu, err := decodeHDU(bytes.NewReader(tt.stream))
			if (err != nil) != tt.wantErr {
				t.Fatalf("decodeHDU() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErrIs != nil && !errors.Is(err, tt.wantErrIs) {
				t.Fatalf("decodeHDU() error = %v, want error matching %v", err, tt.wantErrIs)
			}
			if tt.wantErr {
				return
			}
			if !bytes.Equal(hdu.Data, tt.wantData) {
				t.Errorf("decodeHDU() data = %v, want %v", hdu.Data, tt.wantData)
			}
		})
	}
}

func TestDecodeHeaders(t *testing.T) {
	tests := map[string]struct {
		stream    []byte
		want      Headers
		wantErr   bool
		wantErrIs error
	}{
		"simple header with comments, history and blank cards": {
			stream: hduBytes([]string{
				card("SIMPLE", "T / standard FITS"),
				raw("COMMENT ignored comment"),
				raw("HISTORY ignored history"),
				raw(""),
				card("BITPIX", "-32"),
				card("NAXIS", "0"),
			}, nil),
			want: Headers{"SIMPLE": "T", "BITPIX": "-32", "NAXIS": "0"},
		},
		"quoted values with escaped apostrophes and comments": {
			stream: hduBytes([]string{
				card("XTENSION", "'IMAGE   '  / image extension"),
				card("EXTNAME", "'It''s   '"),
			}, nil),
			want: Headers{"XTENSION": "IMAGE", "EXTNAME": "It's"},
		},
		"duplicate key keeps last value": {
			stream: hduBytes([]string{
				card("SIMPLE", "T"),
				card("NAXIS", "0"),
				card("NAXIS", "2"),
			}, nil),
			want: Headers{"SIMPLE": "T", "NAXIS": "2"},
		},
		"card without equal sign is skipped": {
			stream: hduBytes([]string{
				card("SIMPLE", "T"),
				raw("FOO     bar"),
			}, nil),
			want: Headers{"SIMPLE": "T"},
		},
		"invalid first card": {
			stream: hduBytes([]string{
				raw("FOO     = 1"),
				card("NAXIS", "0"),
			}, nil),
			wantErr: true,
		},
		"empty stream": {
			stream:    nil,
			wantErr:   true,
			wantErrIs: io.EOF,
		},
		"blank first block": {
			stream:    bytes.Repeat([]byte{' '}, BLOCK_SIZE),
			wantErr:   true,
			wantErrIs: io.EOF,
		},
		"header without END card": {
			stream:    padBlock([]byte(card("SIMPLE", "T")+card("NAXIS", "0")), ' '),
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
		"partial block": {
			stream:    []byte(card("SIMPLE", "T")),
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := decodeHeaders(bytes.NewReader(tt.stream))
			if (err != nil) != tt.wantErr {
				t.Fatalf("decodeHeaders() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErrIs != nil && !errors.Is(err, tt.wantErrIs) {
				t.Fatalf("decodeHeaders() error = %v, want error matching %v", err, tt.wantErrIs)
			}
			if tt.wantErr {
				return
			}
			if !maps.Equal(got, tt.want) {
				t.Errorf("decodeHeaders() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestSplitValue(t *testing.T) {
	tests := map[string]struct {
		line string
		want string
	}{
		"value without comment": {
			line: card("NAXIS1", "  10"),
			want: "  10",
		},
		"comment is dropped": {
			line: card("KEY", "42 / comment"),
			want: "42",
		},
		"slash inside quotes is kept": {
			line: card("KEY", "'A/B'  / comment"),
			want: "'A/B'",
		},
		"empty value": {
			line: raw("KEY     ="),
			want: "",
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got := strings.TrimRight(splitValue(tt.line), " ")
			if got != tt.want {
				t.Errorf("splitValue() = %q, want %q", got, tt.want)
			}
		})
	}
}

func TestDecodeValue(t *testing.T) {
	tests := map[string]struct {
		raw  string
		want string
	}{
		"number surrounded by spaces": {
			raw:  "   42  ",
			want: "42",
		},
		"boolean": {
			raw:  "T",
			want: "T",
		},
		"quoted string with trailing spaces": {
			raw:  "'ABC     '",
			want: "ABC",
		},
		"escaped apostrophe": {
			raw:  "'O''Neil'",
			want: "O'Neil",
		},
		"empty quoted string": {
			raw:  "''",
			want: "",
		},
		"quoted spaces only": {
			raw:  "'   '",
			want: "",
		},
		"empty value": {
			raw:  "",
			want: "",
		},
		"unterminated quote": {
			raw:  "'abc",
			want: "abc",
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got := decodeValue(tt.raw)
			if got != tt.want {
				t.Errorf("decodeValue(%q) = %q, want %q", tt.raw, got, tt.want)
			}
		})
	}
}

func TestDecodeBody(t *testing.T) {
	tests := map[string]struct {
		headers   Headers
		stream    []byte
		want      []byte
		wantErr   bool
		wantErrIs error
	}{
		"no data when NAXIS is zero": {
			headers: Headers{"NAXIS": "0"},
			stream:  nil,
			want:    nil,
		},
		"8-bit data padded to block": {
			headers: Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "3"},
			stream:  padBlock([]byte{1, 2, 3}, 0),
			want:    []byte{1, 2, 3},
		},
		"16-bit data swapped to little endian": {
			headers: Headers{"BITPIX": "16", "NAXIS": "1", "NAXIS1": "2"},
			stream:  padBlock([]byte{1, 2, 3, 4}, 0),
			want:    []byte{2, 1, 4, 3},
		},
		"-64-bit data swapped to little endian": {
			headers: Headers{"BITPIX": "-64", "NAXIS": "1", "NAXIS1": "1"},
			stream:  padBlock([]byte{1, 2, 3, 4, 5, 6, 7, 8}, 0),
			want:    []byte{8, 7, 6, 5, 4, 3, 2, 1},
		},
		"data spanning multiple blocks": {
			headers: Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "3000"},
			stream:  padBlock(seqBytes(3000), 0),
			want:    seqBytes(3000),
		},
		"empty stream with required data": {
			headers:   Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "10"},
			stream:    nil,
			wantErr:   true,
			wantErrIs: io.EOF,
		},
		"partial block": {
			headers:   Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "10"},
			stream:    []byte{1, 2},
			wantErr:   true,
			wantErrIs: io.ErrUnexpectedEOF,
		},
		"invalid BITPIX": {
			headers: Headers{"BITPIX": "7", "NAXIS": "1", "NAXIS1": "2"},
			stream:  padBlock([]byte{1, 2}, 0),
			wantErr: true,
		},
		"missing NAXIS": {
			headers: Headers{"BITPIX": "8"},
			stream:  nil,
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := decodeBody(bytes.NewReader(tt.stream), tt.headers)
			if (err != nil) != tt.wantErr {
				t.Fatalf("decodeBody() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErrIs != nil && !errors.Is(err, tt.wantErrIs) {
				t.Fatalf("decodeBody() error = %v, want error matching %v", err, tt.wantErrIs)
			}
			if tt.wantErr {
				return
			}
			if !bytes.Equal(got, tt.want) {
				t.Errorf("decodeBody() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestSwapEndian(t *testing.T) {
	tests := map[string]struct {
		data []byte
		size int
		want []byte
	}{
		"2-byte elements": {
			data: []byte{1, 2, 3, 4},
			size: 2,
			want: []byte{2, 1, 4, 3},
		},
		"4-byte elements": {
			data: []byte{1, 2, 3, 4, 5, 6, 7, 8},
			size: 4,
			want: []byte{4, 3, 2, 1, 8, 7, 6, 5},
		},
		"trailing partial element is untouched": {
			data: []byte{1, 2, 3},
			size: 2,
			want: []byte{2, 1, 3},
		},
		"size one changes nothing": {
			data: []byte{1, 2},
			size: 1,
			want: []byte{1, 2},
		},
		"empty data": {
			data: nil,
			size: 2,
			want: nil,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			swapEndian(tt.data, tt.size)
			if !bytes.Equal(tt.data, tt.want) {
				t.Errorf("swapEndian() = %v, want %v", tt.data, tt.want)
			}
		})
	}
}

func TestIsBlank(t *testing.T) {
	tests := map[string]struct {
		b    []byte
		want bool
	}{
		"empty slice": {
			b:    nil,
			want: true,
		},
		"zeros and spaces": {
			b:    []byte{0, ' ', 0},
			want: true,
		},
		"contains letter": {
			b:    []byte("  A "),
			want: false,
		},
		"contains tab": {
			b:    []byte{' ', '\t'},
			want: false,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			if got := isBlank(tt.b); got != tt.want {
				t.Errorf("isBlank() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestIntHeader(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		key     string
		want    int
		wantErr bool
	}{
		"valid positive value": {
			headers: Headers{"NAXIS": "2"},
			key:     "NAXIS",
			want:    2,
		},
		"valid negative value": {
			headers: Headers{"BITPIX": "-32"},
			key:     "BITPIX",
			want:    -32,
		},
		"missing key": {
			headers: Headers{},
			key:     "NAXIS",
			wantErr: true,
		},
		"not a number": {
			headers: Headers{"NAXIS": "abc"},
			key:     "NAXIS",
			wantErr: true,
		},
		"empty value": {
			headers: Headers{"NAXIS": ""},
			key:     "NAXIS",
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := intHeader(tt.headers, tt.key)
			if (err != nil) != tt.wantErr {
				t.Fatalf("intHeader() error = %v, wantErr %v", err, tt.wantErr)
			}
			if got != tt.want {
				t.Errorf("intHeader() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestDataSize(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    int
		wantErr bool
	}{
		"missing NAXIS": {
			headers: Headers{"BITPIX": "8"},
			wantErr: true,
		},
		"invalid NAXIS": {
			headers: Headers{"NAXIS": "x"},
			wantErr: true,
		},
		"negative NAXIS": {
			headers: Headers{"NAXIS": "-1"},
			wantErr: true,
		},
		"NAXIS above maximum": {
			headers: Headers{"NAXIS": "1000"},
			wantErr: true,
		},
		"NAXIS zero": {
			headers: Headers{"NAXIS": "0"},
			want:    0,
		},
		"missing BITPIX": {
			headers: Headers{"NAXIS": "1", "NAXIS1": "10"},
			wantErr: true,
		},
		"invalid BITPIX": {
			headers: Headers{"BITPIX": "12", "NAXIS": "1", "NAXIS1": "10"},
			wantErr: true,
		},
		"missing NAXISn": {
			headers: Headers{"BITPIX": "8", "NAXIS": "2", "NAXIS1": "3"},
			wantErr: true,
		},
		"NAXISn not positive": {
			headers: Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "0"},
			wantErr: true,
		},
		"NAXISn not a number": {
			headers: Headers{"BITPIX": "8", "NAXIS": "1", "NAXIS1": "abc"},
			wantErr: true,
		},
		"8-bit 2D image": {
			headers: Headers{"BITPIX": "8", "NAXIS": "2", "NAXIS1": "3", "NAXIS2": "4"},
			want:    12,
		},
		"16-bit 1D data": {
			headers: Headers{"BITPIX": "16", "NAXIS": "1", "NAXIS1": "5"},
			want:    10,
		},
		"-64-bit 2D image": {
			headers: Headers{"BITPIX": "-64", "NAXIS": "2", "NAXIS1": "2", "NAXIS2": "3"},
			want:    48,
		},
		"dimension product overflows limit": {
			headers: Headers{"BITPIX": "8", "NAXIS": "2", "NAXIS1": "2147483647", "NAXIS2": "2147483647"},
			wantErr: true,
		},
		"size exceeds limit after BITPIX scaling": {
			headers: Headers{"BITPIX": "-64", "NAXIS": "1", "NAXIS1": "576460752303423488"},
			wantErr: true,
		},
		"size exactly at limit": {
			headers: Headers{"BITPIX": "-64", "NAXIS": "1", "NAXIS1": "144115188075855872"},
			want:    1 << 60,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := dataSize(tt.headers)
			if (err != nil) != tt.wantErr {
				t.Fatalf("dataSize() error = %v, wantErr %v", err, tt.wantErr)
			}
			if got != tt.want {
				t.Errorf("dataSize() = %d, want %d", got, tt.want)
			}
		})
	}
}

// raw returns s padded with spaces to the length of a single FITS header card.
func raw(s string) string {
	return fmt.Sprintf("%-*s", HDU_LINE_LENGTH, s)
}

// card returns a header card with the given keyword and value.
func card(key, value string) string {
	return raw(fmt.Sprintf("%-8s= %s", key, value))
}

// hduBytes builds a complete HDU: header cards followed by an END card, padded
// with spaces to BLOCK_SIZE, and then data padded with zeros to BLOCK_SIZE.
func hduBytes(cards []string, data []byte) []byte {
	lines := append(append([]string{}, cards...), raw("END"))
	out := padBlock([]byte(strings.Join(lines, "")), ' ')
	if len(data) == 0 {
		return out
	}
	out = append(out, data...)
	return padBlock(out, 0)
}

// padBlock pads b with fill bytes up to a multiple of BLOCK_SIZE.
func padBlock(b []byte, fill byte) []byte {
	rem := len(b) % BLOCK_SIZE
	if rem == 0 {
		return b
	}
	return append(b, bytes.Repeat([]byte{fill}, BLOCK_SIZE-rem)...)
}

// seqBytes returns n bytes with values 0, 1, 2, ... wrapping at 255.
func seqBytes(n int) []byte {
	b := make([]byte, n)
	for i := range b {
		b[i] = byte(i)
	}
	return b
}
