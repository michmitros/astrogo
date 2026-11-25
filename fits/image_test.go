package fits

import (
	"encoding/binary"
	"image"
	"maps"
	"math"
	"slices"
	"strconv"
	"testing"
)

func TestToImage(t *testing.T) {
	tests := map[string]struct {
		hdu      HDU
		contrast float64
		samples  int
		want     []uint16
		wantErr  bool
	}{
		"8-bit image flipped vertically": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"16-bit image": {
			hdu:      newImageHDU(16, 2, 2, encodePixels(16, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"32-bit image": {
			hdu:      newImageHDU(32, 2, 2, encodePixels(32, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"-32-bit image": {
			hdu:      newImageHDU(-32, 2, 2, encodePixels(-32, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"-64-bit image": {
			hdu:      newImageHDU(-64, 2, 2, encodePixels(-64, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"BZERO and negative BSCALE applied": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{0, 0, 1, 1}), Headers{"BZERO": "0", "BSCALE": "-1"}),
			contrast: 1,
			want:     []uint16{0, 0, 65535, 65535},
		},
		"invalid BSCALE falls back to default": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{0, 0, 255, 255}), Headers{"BSCALE": "abc"}),
			contrast: 1,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"constant image renders black": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{7, 7, 7, 7}), nil),
			contrast: 1,
			want:     []uint16{0, 0, 0, 0},
		},
		"NaN and infinite pixels render black": {
			hdu:      newImageHDU(-32, 2, 2, encodePixels(-32, []float64{math.NaN(), math.Inf(1), 0, 255}), nil),
			contrast: 1,
			want:     []uint16{0, 65535, 0, 0},
		},
		"sampling limit is applied": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{0, 0, 255, 255}), nil),
			contrast: 1,
			samples:  2,
			want:     []uint16{65535, 65535, 0, 0},
		},
		"unsupported BITPIX reads zeros": {
			hdu:      newImageHDU(7, 2, 2, []byte{1, 2, 3, 4}, nil),
			contrast: 1,
			want:     []uint16{0, 0, 0, 0},
		},
		"NAXIS is not 2": {
			hdu:      newImageHDU(8, 2, 2, encodePixels(8, []float64{0, 0, 255, 255}), Headers{"NAXIS": "3"}),
			contrast: 1,
			wantErr:  true,
		},
		"NAXIS is missing": {
			hdu:      HDU{Headers: Headers{"BITPIX": "8", "NAXIS1": "2", "NAXIS2": "2"}, Data: []byte{0, 0, 255, 255}},
			contrast: 1,
			wantErr:  true,
		},
		"NAXIS1 is zero": {
			hdu:      newImageHDU(8, 0, 2, nil, nil),
			contrast: 1,
			wantErr:  true,
		},
		"NAXIS2 is negative": {
			hdu:      newImageHDU(8, 2, -1, nil, nil),
			contrast: 1,
			wantErr:  true,
		},
		"data too short": {
			hdu:      newImageHDU(8, 2, 2, []byte{1, 2, 3}, nil),
			contrast: 1,
			wantErr:  true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			img, err := ToImage(tt.hdu, tt.contrast, tt.samples)
			if (err != nil) != tt.wantErr {
				t.Fatalf("ToImage() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErr {
				return
			}
			if got := pixelsOf(img); !slices.Equal(got, tt.want) {
				t.Errorf("ToImage() pixels = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestZscale(t *testing.T) {
	tests := map[string]struct {
		values   []float64
		contrast float64
		samples  int
		wantLo   float64
		wantHi   float64
	}{
		"empty input": {
			values:   nil,
			contrast: 1,
			wantLo:   0,
			wantHi:   0,
		},
		"only NaN and infinite values": {
			values:   []float64{math.NaN(), math.Inf(-1), math.Inf(1)},
			contrast: 1,
			wantLo:   0,
			wantHi:   0,
		},
		"too few samples returns min and max": {
			values:   []float64{1, 2, 3, 4},
			contrast: 1,
			wantLo:   1,
			wantHi:   4,
		},
		"constant values give empty range": {
			values:   slices.Repeat([]float64{5}, 100),
			contrast: 1,
			wantLo:   5,
			wantHi:   5,
		},
		"linear ramp gives narrowed range": {
			values:   linspace(0, 999, 1000),
			contrast: 1,
			wantLo:   0,
			wantHi:   999,
		},
		"larger contrast gives narrower range": {
			values:   linspace(0, 999, 1000),
			contrast: 10,
			wantLo:   450,
			wantHi:   549.9,
		},
		"outliers are rejected by iterative fit": {
			values:   append(linspace(0, 99, 100), 1e9, -1e9),
			contrast: 1,
			wantLo:   -1,
			wantHi:   100,
		},
		"sampling reduces number of pixels": {
			values:   linspace(0, 999, 1000),
			contrast: 1,
			samples:  100,
			wantLo:   0,
			wantHi:   990,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			lo, hi := zscale(slices.Clone(tt.values), tt.contrast, tt.samples)
			if !closeTo(lo, tt.wantLo) || !closeTo(hi, tt.wantHi) {
				t.Errorf("zscale() = (%v, %v), want (%v, %v)", lo, hi, tt.wantLo, tt.wantHi)
			}
		})
	}
}

func TestReadValue(t *testing.T) {
	tests := map[string]struct {
		bitpix int
		data   []byte
		i      int
		want   float64
	}{
		"8-bit unsigned": {
			bitpix: 8,
			data:   []byte{200},
			i:      0,
			want:   200,
		},
		"16-bit signed negative": {
			bitpix: 16,
			data:   encodePixels(16, []float64{-2}),
			i:      0,
			want:   -2,
		},
		"32-bit signed second element": {
			bitpix: 32,
			data:   encodePixels(32, []float64{1, -70000}),
			i:      1,
			want:   -70000,
		},
		"-32-bit float": {
			bitpix: -32,
			data:   encodePixels(-32, []float64{1.5}),
			i:      0,
			want:   1.5,
		},
		"-64-bit float": {
			bitpix: -64,
			data:   encodePixels(-64, []float64{-2.25}),
			i:      0,
			want:   -2.25,
		},
		"unsupported BITPIX returns zero": {
			bitpix: 4,
			data:   []byte{1, 2, 3, 4},
			i:      0,
			want:   0,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			if got := readValue(tt.data, tt.bitpix, tt.i); got != tt.want {
				t.Errorf("readValue() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestFloatHeader(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		key     string
		def     float64
		want    float64
	}{
		"valid value": {
			headers: Headers{"BZERO": "32768.5"},
			key:     "BZERO",
			def:     0,
			want:    32768.5,
		},
		"missing key returns default": {
			headers: Headers{},
			key:     "BSCALE",
			def:     1,
			want:    1,
		},
		"invalid value returns default": {
			headers: Headers{"BSCALE": "abc"},
			key:     "BSCALE",
			def:     1,
			want:    1,
		},
		"empty value returns default": {
			headers: Headers{"BSCALE": ""},
			key:     "BSCALE",
			def:     1,
			want:    1,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			if got := floatHeader(tt.headers, tt.key, tt.def); got != tt.want {
				t.Errorf("floatHeader() = %v, want %v", got, tt.want)
			}
		})
	}
}

// newImageHDU returns a 2D image HDU with the given BITPIX, dimensions, data and
// extra header cards that override the generated ones.
func newImageHDU(bitpix, w, h int, data []byte, extra Headers) HDU {
	headers := Headers{
		"XTENSION": string(ExtImage),
		"BITPIX":   strconv.Itoa(bitpix),
		"NAXIS":    "2",
		"NAXIS1":   strconv.Itoa(w),
		"NAXIS2":   strconv.Itoa(h),
	}
	maps.Copy(headers, extra)
	return HDU{Headers: headers, Data: data}
}

// encodePixels encodes values as little endian raw data of the given BITPIX.
func encodePixels(bitpix int, values []float64) []byte {
	var out []byte
	for _, v := range values {
		switch bitpix {
		case 8:
			out = append(out, byte(v))
		case 16:
			out = binary.LittleEndian.AppendUint16(out, uint16(int16(v)))
		case 32:
			out = binary.LittleEndian.AppendUint32(out, uint32(int32(v)))
		case -32:
			out = binary.LittleEndian.AppendUint32(out, math.Float32bits(float32(v)))
		case -64:
			out = binary.LittleEndian.AppendUint64(out, math.Float64bits(v))
		}
	}
	return out
}

// pixelsOf returns the gray values of img in row-major order, starting at the top row.
func pixelsOf(img *image.Gray16) []uint16 {
	b := img.Bounds()
	var out []uint16
	for y := b.Min.Y; y < b.Max.Y; y++ {
		for x := b.Min.X; x < b.Max.X; x++ {
			out = append(out, img.Gray16At(x, y).Y)
		}
	}
	return out
}

// linspace returns n evenly spaced values from start to end inclusive.
func linspace(start, end float64, n int) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = start + (end-start)*float64(i)/float64(n-1)
	}
	return out
}

// closeTo reports whether a and b differ by less than a small tolerance.
func closeTo(a, b float64) bool {
	return math.Abs(a-b) < 1e-6
}
