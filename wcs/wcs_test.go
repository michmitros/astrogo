package wcs

import (
	"math"
	"testing"

	"github.com/michmitros/astrogo/fits"
)

const (
	matrixTol = 1e-12
	skyTol    = 1e-5 // deg
)

func TestNewFromHeaders(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		wantCD  [2][2]float64
		wantRef [2]float64
		wantVal [2]float64
		wantErr bool
	}{
		"CDELT with default PC matrix": {
			headers: newHeaders(nil),
			wantCD:  [2][2]float64{{-0.001, 0}, {0, 0.001}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"PC matrix scales CDELT": {
			headers: newHeaders(fits.Headers{"PC1_1": "0.6", "PC1_2": "-0.8", "PC2_1": "0.8", "PC2_2": "0.6"}),
			wantCD:  [2][2]float64{{-0.0006, 0.0008}, {0.0008, 0.0006}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"CD matrix overrides CDELT": {
			headers: newHeaders(fits.Headers{"CD1_1": "-0.0006", "CD1_2": "0.0008", "CD2_1": "0.0008", "CD2_2": "0.0006"}),
			wantCD:  [2][2]float64{{-0.0006, 0.0008}, {0.0008, 0.0006}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"missing CD elements are zero": {
			headers: newHeaders(fits.Headers{"CD1_1": "-0.001", "CD2_2": "0.001"}),
			wantCD:  [2][2]float64{{-0.001, 0}, {0, 0.001}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"D exponent in CDELT": {
			headers: newHeaders(fits.Headers{"CDELT1": "-1.0D-03", "CDELT2": "1.0d-03"}),
			wantCD:  [2][2]float64{{-0.001, 0}, {0, 0.001}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"CDELT defaults to one": {
			headers: newHeaders(nil, "CDELT1", "CDELT2"),
			wantCD:  [2][2]float64{{1, 0}, {0, 1}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"CRPIX and CRVAL default to zero": {
			headers: newHeaders(nil, "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2"),
			wantCD:  [2][2]float64{{-0.001, 0}, {0, 0.001}},
			wantRef: [2]float64{0, 0},
			wantVal: [2]float64{0, 0},
		},
		"CUNIT deg accepted case insensitive with spaces": {
			headers: newHeaders(fits.Headers{"CUNIT1": "deg", "CUNIT2": " DEG "}),
			wantCD:  [2][2]float64{{-0.001, 0}, {0, 0.001}},
			wantRef: [2]float64{50.5, 50.5},
			wantVal: [2]float64{150, 2},
		},
		"missing NAXIS1": {
			headers: newHeaders(nil, "NAXIS1"),
			wantErr: true,
		},
		"missing NAXIS2": {
			headers: newHeaders(nil, "NAXIS2"),
			wantErr: true,
		},
		"invalid NAXIS1": {
			headers: newHeaders(fits.Headers{"NAXIS1": "abc"}),
			wantErr: true,
		},
		"zero NAXIS1": {
			headers: newHeaders(fits.Headers{"NAXIS1": "0"}),
			wantErr: true,
		},
		"negative NAXIS2": {
			headers: newHeaders(fits.Headers{"NAXIS2": "-5"}),
			wantErr: true,
		},
		"invalid CRPIX1": {
			headers: newHeaders(fits.Headers{"CRPIX1": "x"}),
			wantErr: true,
		},
		"NaN CRVAL2": {
			headers: newHeaders(fits.Headers{"CRVAL2": "NaN"}),
			wantErr: true,
		},
		"missing CTYPE1": {
			headers: newHeaders(nil, "CTYPE1"),
			wantErr: true,
		},
		"unsupported projection SIN": {
			headers: newHeaders(fits.Headers{"CTYPE2": "DEC--SIN"}),
			wantErr: true,
		},
		"CTYPE too short": {
			headers: newHeaders(fits.Headers{"CTYPE1": "RA-TAN"}),
			wantErr: true,
		},
		"CTYPE without dash at position 5": {
			headers: newHeaders(fits.Headers{"CTYPE1": "RAXXXTAN"}),
			wantErr: true,
		},
		"unsupported unit": {
			headers: newHeaders(fits.Headers{"CUNIT1": "rad"}),
			wantErr: true,
		},
		"unsupported unit on second axis": {
			headers: newHeaders(fits.Headers{"CUNIT2": "arcsec"}),
			wantErr: true,
		},
		"invalid CDELT1": {
			headers: newHeaders(fits.Headers{"CDELT1": "abc"}),
			wantErr: true,
		},
		"NaN PC2_1": {
			headers: newHeaders(fits.Headers{"PC2_1": "NaN"}),
			wantErr: true,
		},
		"out of range CD element": {
			headers: newHeaders(fits.Headers{"CD1_1": "1e400"}),
			wantErr: true,
		},
		"singular matrix from CDELT": {
			headers: newHeaders(fits.Headers{"CDELT1": "0"}),
			wantErr: true,
		},
		"singular matrix from CD": {
			headers: newHeaders(fits.Headers{"CD1_1": "1", "CD1_2": "2", "CD2_1": "2", "CD2_2": "4"}),
			wantErr: true,
		},
		"determinant overflows to infinity": {
			headers: newHeaders(fits.Headers{"CD1_1": "1e200", "CD1_2": "0", "CD2_1": "0", "CD2_2": "1e200"}),
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w, err := NewFromHeaders(tt.headers)

			if (err != nil) != tt.wantErr {
				t.Fatalf("NewFromHeaders() error = %v, wantErr %v", err, tt.wantErr)
			}
			if tt.wantErr {
				return
			}

			if !approxMatrix(w.cd, tt.wantCD) {
				t.Errorf("NewFromHeaders() cd = %v, want %v", w.cd, tt.wantCD)
			}

			if w.crpix != tt.wantRef {
				t.Errorf("NewFromHeaders() crpix = %v, want %v", w.crpix, tt.wantRef)
			}

			if w.crval != tt.wantVal {
				t.Errorf("NewFromHeaders() crval = %v, want %v", w.crval, tt.wantVal)
			}

			if w.naxis != [2]int{100, 100} {
				t.Errorf("NewFromHeaders() naxis = %v, want [100 100]", w.naxis)
			}
		})
	}
}

func TestInverseMatrix(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
	}{
		"diagonal matrix": {
			headers: newHeaders(nil),
		},
		"rotated matrix": {
			headers: newHeaders(fits.Headers{"PC1_1": "0.6", "PC1_2": "-0.8", "PC2_1": "0.8", "PC2_2": "0.6"}),
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, tt.headers)

			// cd * cdInv must be the identity matrix.
			var prod [2][2]float64
			for i := range 2 {
				for j := range 2 {
					for k := range 2 {
						prod[i][j] += w.cd[i][k] * w.cdInv[k][j]
					}
				}
			}

			want := [2][2]float64{{1, 0}, {0, 1}}
			if !approxMatrix(prod, want) {
				t.Errorf("cd * cdInv = %v, want identity", prod)
			}
		})
	}
}

func TestPixelToWorld(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		pixel   Pixel
		want    SkyCoord
	}{
		"reference pixel": {
			headers: newHeaders(fits.Headers{"CRPIX1": "1", "CRPIX2": "1"}),
			pixel:   Pixel{0, 0},
			want:    SkyCoord{RA: 150, Dec: 2},
		},
		"pixel left and below of center": {
			headers: newHeaders(nil),
			pixel:   Pixel{49, 49},
			want:    SkyCoord{RA: 150.0005, Dec: 1.9995},
		},
		"pixel right and above of center": {
			headers: newHeaders(nil),
			pixel:   Pixel{50, 50},
			want:    SkyCoord{RA: 149.9995, Dec: 2.0005},
		},
		"RA wraps into [0, 360)": {
			headers: newHeaders(fits.Headers{"CRVAL1": "0", "CRPIX1": "1", "CRPIX2": "1"}),
			pixel:   Pixel{1, 0},
			want:    SkyCoord{RA: 359.999, Dec: 2},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, tt.headers)

			got, err := w.PixelToWorld(tt.pixel)
			if err != nil {
				t.Fatalf("PixelToWorld() error = %v", err)
			}

			if !approxSkyCoord(*got, tt.want, skyTol) {
				t.Errorf("PixelToWorld(%v) = %v, want %v", tt.pixel, *got, tt.want)
			}
			if got.RA < 0 || got.RA >= 360 {
				t.Errorf("PixelToWorld(%v) RA = %v, want within [0, 360)", tt.pixel, got.RA)
			}
		})
	}
}

func TestPixelsToWorld(t *testing.T) {
	tests := map[string]struct {
		pixels []Pixel
		want   []SkyCoord
	}{
		"no pixels": {
			pixels: nil,
			want:   []SkyCoord{},
		},
		"multiple pixels keep order": {
			pixels: []Pixel{{49, 49}, {50, 50}},
			want:   []SkyCoord{{RA: 150.0005, Dec: 1.9995}, {RA: 149.9995, Dec: 2.0005}},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))

			got, err := w.PixelsToWorld(tt.pixels)
			if err != nil {
				t.Fatalf("PixelsToWorld() error = %v", err)
			}

			if len(got) != len(tt.want) {
				t.Fatalf("PixelsToWorld() returned %d coords, want %d", len(got), len(tt.want))
			}
			for i := range got {
				if !approxSkyCoord(got[i], tt.want[i], skyTol) {
					t.Errorf("PixelsToWorld() #%d = %v, want %v", i, got[i], tt.want[i])
				}
			}
		})
	}
}

func TestWorldToPixel(t *testing.T) {
	tests := map[string]struct {
		x, y float64
		want Pixel
	}{
		"integer pixel center": {
			x:    49,
			y:    49,
			want: Pixel{49, 49},
		},
		"fraction below half rounds down": {
			x:    10.4,
			y:    20.4,
			want: Pixel{10, 20},
		},
		"fraction above half rounds up": {
			x:    10.6,
			y:    20.6,
			want: Pixel{11, 21},
		},
		"negative fraction rounds to -1": {
			x:    -0.6,
			y:    0,
			want: Pixel{-1, 0},
		},
		"negative fraction rounds to 0": {
			x:    -0.4,
			y:    0,
			want: Pixel{0, 0},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))

			coord := w.pixelToWorld(tt.x, tt.y)
			got, err := w.WorldToPixel(coord)
			if err != nil {
				t.Fatalf("WorldToPixel() error = %v", err)
			}

			if got != tt.want {
				t.Errorf("WorldToPixel(%v) = %v, want %v", coord, got, tt.want)
			}
		})
	}
}

func TestWorldToPixelErrors(t *testing.T) {
	tests := map[string]struct {
		coord SkyCoord
	}{
		"NaN RA": {
			coord: SkyCoord{RA: math.NaN(), Dec: 2},
		},
		"infinite RA": {
			coord: SkyCoord{RA: math.Inf(1), Dec: 2},
		},
		"declination above 90": {
			coord: SkyCoord{RA: 150, Dec: 91},
		},
		"NaN declination": {
			coord: SkyCoord{RA: 150, Dec: math.NaN()},
		},
		"point beyond 90 degrees from reference": {
			coord: SkyCoord{RA: 330, Dec: -2},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))

			if _, err := w.WorldToPixel(tt.coord); err == nil {
				t.Errorf("WorldToPixel(%v) expected an error, got none", tt.coord)
			}

			if _, _, err := w.worldToPixel(tt.coord); err == nil {
				t.Errorf("worldToPixel(%v) expected an error, got none", tt.coord)
			}

			if w.Contains(tt.coord) {
				t.Errorf("Contains(%v) = true, want false for unprojectable coordinate", tt.coord)
			}
		})
	}
}

func TestWorldToPixels(t *testing.T) {
	tests := map[string]struct {
		pixels  []Pixel
		invalid bool
	}{
		"no coordinates": {
			pixels: nil,
		},
		"all coordinates valid": {
			pixels: []Pixel{{0, 0}, {49, 49}, {99, 99}},
		},
		"invalid coordinate at the end": {
			pixels:  []Pixel{{0, 0}, {49, 49}},
			invalid: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))

			coords, err := w.PixelsToWorld(tt.pixels)
			if err != nil {
				t.Fatalf("PixelsToWorld() error = %v", err)
			}

			if tt.invalid {
				coords = append(coords, SkyCoord{RA: 330, Dec: -2})
			}

			got, err := w.WorldToPixels(coords)
			if (err != nil) != tt.invalid {
				t.Fatalf("WorldToPixels() error = %v, wantErr %v", err, tt.invalid)
			}

			if tt.invalid {
				if got != nil {
					t.Errorf("WorldToPixels() = %v, want nil on error", got)
				}

				return
			}

			if len(got) != len(tt.pixels) {
				t.Fatalf("WorldToPixels() returned %d pixels, want %d", len(got), len(tt.pixels))
			}
			for i := range got {
				if got[i] != tt.pixels[i] {
					t.Errorf("WorldToPixels() #%d = %v, want %v", i, got[i], tt.pixels[i])
				}
			}
		})
	}
}

func TestWorldToPixelRoundTrip(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		pixels  []Pixel
	}{
		"CDELT with no rotation": {
			headers: newHeaders(nil),
			pixels:  []Pixel{{0, 0}, {99, 99}, {10, 90}, {49, 49}, {-5, 120}},
		},
		"rotated PC matrix": {
			headers: newHeaders(fits.Headers{
				"PC1_1": "0.8660254037844387", "PC1_2": "-0.5",
				"PC2_1": "0.5", "PC2_2": "0.8660254037844387",
			}),
			pixels: []Pixel{{0, 0}, {99, 99}, {10, 90}, {49, 49}, {-5, 120}},
		},
		"CD matrix with skew": {
			headers: newHeaders(fits.Headers{
				"CD1_1": "-0.001", "CD1_2": "0.0002",
				"CD2_1": "0.0001", "CD2_2": "0.0011",
				"CRPIX1": "20", "CRPIX2": "70",
			}),
			pixels: []Pixel{{0, 0}, {19, 69}, {60, 10}, {99, 99}},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, tt.headers)

			for _, p := range tt.pixels {
				coord := w.pixelToWorld(float64(p[0]), float64(p[1]))

				x, y, err := w.worldToPixel(coord)
				if err != nil {
					t.Fatalf("worldToPixel(%v) error = %v", coord, err)
				}

				if !approxEqual(x, float64(p[0]), 1e-6) || !approxEqual(y, float64(p[1]), 1e-6) {
					t.Errorf("worldToPixel(pixelToWorld(%v)) = (%g, %g), want %v", p, x, y, p)
				}

				got, err := w.WorldToPixel(coord)
				if err != nil {
					t.Fatalf("WorldToPixel(%v) error = %v", coord, err)
				}

				if got != p {
					t.Errorf("WorldToPixel(pixelToWorld(%v)) = %v, want %v", p, got, p)
				}
			}
		})
	}
}

func TestFootprint(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		corners [4][2]float64
	}{
		"square image": {
			headers: newHeaders(nil),
			corners: [4][2]float64{{-0.5, -0.5}, {99.5, -0.5}, {99.5, 99.5}, {-0.5, 99.5}},
		},
		"rectangular image": {
			headers: newHeaders(fits.Headers{"NAXIS1": "200", "NAXIS2": "50"}),
			corners: [4][2]float64{{-0.5, -0.5}, {199.5, -0.5}, {199.5, 49.5}, {-0.5, 49.5}},
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, tt.headers)
			got := w.Footprint()

			if len(got) != len(tt.corners) {
				t.Fatalf("Footprint() returned %d corners, want %d", len(got), len(tt.corners))
			}
			for i, c := range tt.corners {
				want := w.pixelToWorld(c[0], c[1])

				if !approxSkyCoord(got[i], want, matrixTol) {
					t.Errorf("Footprint() corner #%d = %v, want %v", i, got[i], want)
				}
			}
		})
	}
}

func TestContains(t *testing.T) {
	tests := map[string]struct {
		x, y float64
		want bool
	}{
		"center of image": {
			x:    49.5,
			y:    49.5,
			want: true,
		},
		"first pixel": {
			x:    0,
			y:    0,
			want: true,
		},
		"last pixel": {
			x:    99,
			y:    99,
			want: true,
		},
		"left of first column": {
			x:    -1,
			y:    0,
			want: false,
		},
		"right of last column": {
			x:    100,
			y:    0,
			want: false,
		},
		"above last row": {
			x:    0,
			y:    100,
			want: false,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))

			coord := w.pixelToWorld(tt.x, tt.y)
			if got := w.Contains(coord); got != tt.want {
				t.Errorf("Contains(%v) = %v, want %v", coord, got, tt.want)
			}
		})
	}
}

func TestPixelScale(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		want    float64
	}{
		"CDELT 0.001 deg per pixel": {
			headers: newHeaders(nil),
			want:    3.6,
		},
		"rotated CD matrix with same determinant": {
			headers: newHeaders(fits.Headers{"CD1_1": "0.0006", "CD1_2": "-0.0008", "CD2_1": "0.0008", "CD2_2": "0.0006"}),
			want:    3.6,
		},
		"CDELT 0.0005 deg per pixel": {
			headers: newHeaders(fits.Headers{"CDELT1": "-0.0005", "CDELT2": "0.0005"}),
			want:    1.8,
		},
		"default CDELT of one degree": {
			headers: newHeaders(nil, "CDELT1", "CDELT2"),
			want:    3600,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, tt.headers)

			got, err := w.PixelScale()
			if err != nil {
				t.Fatalf("PixelScale() error = %v", err)
			}

			if !approxEqual(got, tt.want, 1e-9) {
				t.Errorf("PixelScale() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestSlice(t *testing.T) {
	tests := map[string]struct {
		xMin, yMin, xMax, yMax int
		wantErr                bool
	}{
		"full image": {
			xMin: 0, yMin: 0, xMax: 100, yMax: 100,
		},
		"sub region": {
			xMin: 10, yMin: 5, xMax: 20, yMax: 15,
		},
		"negative xMin": {
			xMin: -1, yMin: 0, xMax: 10, yMax: 10,
			wantErr: true,
		},
		"negative yMin": {
			xMin: 0, yMin: -1, xMax: 10, yMax: 10,
			wantErr: true,
		},
		"empty column range": {
			xMin: 10, yMin: 0, xMax: 10, yMax: 10,
			wantErr: true,
		},
		"empty row range": {
			xMin: 0, yMin: 10, xMax: 10, yMax: 5,
			wantErr: true,
		},
		"xMax beyond image": {
			xMin: 0, yMin: 0, xMax: 101, yMax: 10,
			wantErr: true,
		},
		"yMax beyond image": {
			xMin: 0, yMin: 0, xMax: 10, yMax: 101,
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			w := newTestWCS(t, newHeaders(nil))
			orig := *w

			s, err := w.Slice(tt.xMin, tt.yMin, tt.xMax, tt.yMax)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Slice() error = %v, wantErr %v", err, tt.wantErr)
			}

			if tt.wantErr {
				if s != nil {
					t.Errorf("Slice() = %v, want nil on error", s)
				}
				return
			}

			wantNaxis := [2]int{tt.xMax - tt.xMin, tt.yMax - tt.yMin}
			if s.naxis != wantNaxis {
				t.Errorf("Slice() naxis = %v, want %v", s.naxis, wantNaxis)
			}

			got := s.pixelToWorld(0, 0)
			want := w.pixelToWorld(float64(tt.xMin), float64(tt.yMin))

			if !approxSkyCoord(got, want, matrixTol) {
				t.Errorf("Slice() pixel (0, 0) = %v, want %v", got, want)
			}

			if *w != orig {
				t.Errorf("Slice() modified the original WCS")
			}
		})
	}
}

func TestHeaderFloat(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		def     float64
		want    float64
		wantErr bool
	}{
		"absent key returns default": {
			headers: fits.Headers{},
			def:     7,
			want:    7,
		},
		"plain decimal with spaces": {
			headers: fits.Headers{"KEY": " 3.5 "},
			want:    3.5,
		},
		"upper case D exponent": {
			headers: fits.Headers{"KEY": "1.0D-03"},
			want:    0.001,
		},
		"lower case d exponent": {
			headers: fits.Headers{"KEY": "2d1"},
			want:    20,
		},
		"invalid text": {
			headers: fits.Headers{"KEY": "abc"},
			wantErr: true,
		},
		"NaN is rejected": {
			headers: fits.Headers{"KEY": "NaN"},
			wantErr: true,
		},
		"infinity is rejected": {
			headers: fits.Headers{"KEY": "Inf"},
			wantErr: true,
		},
		"out of range value is rejected": {
			headers: fits.Headers{"KEY": "1e400"},
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := headerFloat(tt.headers, "KEY", tt.def)
			if (err != nil) != tt.wantErr {
				t.Fatalf("headerFloat() error = %v, wantErr %v", err, tt.wantErr)
			}

			if got != tt.want {
				t.Errorf("headerFloat() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestHeaderInt(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		want    int
		wantErr bool
	}{
		"valid negative value": {
			headers: fits.Headers{"KEY": "-32"},
			want:    -32,
		},
		"value with spaces": {
			headers: fits.Headers{"KEY": " 7 "},
			want:    7,
		},
		"missing key": {
			headers: fits.Headers{},
			wantErr: true,
		},
		"empty value": {
			headers: fits.Headers{"KEY": ""},
			wantErr: true,
		},
		"non integer value": {
			headers: fits.Headers{"KEY": "1.5"},
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := headerInt(tt.headers, "KEY")
			if (err != nil) != tt.wantErr {
				t.Fatalf("headerInt() error = %v, wantErr %v", err, tt.wantErr)
			}

			if got != tt.want {
				t.Errorf("headerInt() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestHeaderString(t *testing.T) {
	tests := map[string]struct {
		headers fits.Headers
		want    string
		wantErr bool
	}{
		"value is trimmed": {
			headers: fits.Headers{"KEY": "  RA---TAN  "},
			want:    "RA---TAN",
		},
		"missing key": {
			headers: fits.Headers{},
			wantErr: true,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			got, err := headerString(tt.headers, "KEY")
			if (err != nil) != tt.wantErr {
				t.Fatalf("headerString() error = %v, wantErr %v", err, tt.wantErr)
			}

			if got != tt.want {
				t.Errorf("headerString() = %q, want %q", got, tt.want)
			}
		})
	}
}

// newHeaders returns a valid TAN header for a 100x100 image centered at (150, 2) degrees,
// with CDELT of -0.001 and 0.001 degrees per pixel. Entries in overrides replace the defaults,
// and keys listed in remove are deleted afterwards.
func newHeaders(overrides fits.Headers, remove ...string) fits.Headers {
	h := fits.Headers{
		"NAXIS1": "100",
		"NAXIS2": "100",
		"CTYPE1": "RA---TAN",
		"CTYPE2": "DEC--TAN",
		"CRPIX1": "50.5",
		"CRPIX2": "50.5",
		"CRVAL1": "150",
		"CRVAL2": "2",
		"CDELT1": "-0.001",
		"CDELT2": "0.001",
	}

	for k, v := range overrides {
		h[k] = v
	}

	for _, k := range remove {
		delete(h, k)
	}

	return h
}

// newTestWCS creates a WCS from headers and fails the test if they are invalid.
func newTestWCS(t *testing.T, h fits.Headers) *WCS {
	t.Helper()

	w, err := NewFromHeaders(h)
	if err != nil {
		t.Fatalf("NewFromHeaders() error = %v", err)
	}

	return w
}

// approxEqual reports whether a and b differ by no more than tol.
func approxEqual(a, b, tol float64) bool {
	return math.Abs(a-b) <= tol
}

// approxSkyCoord reports whether both coordinates are equal within tol degrees.
func approxSkyCoord(a, b SkyCoord, tol float64) bool {
	return approxEqual(a.RA, b.RA, tol) && approxEqual(a.Dec, b.Dec, tol)
}

// approxMatrix reports whether two 2x2 matrices are equal within matrixTol.
func approxMatrix(a, b [2][2]float64) bool {
	for i := range 2 {
		for j := range 2 {
			if !approxEqual(a[i][j], b[i][j], matrixTol) {
				return false
			}
		}
	}

	return true
}
