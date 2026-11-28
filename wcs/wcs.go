package wcs

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/michmitros/astrogo/fits"
)

const (
	deg2rad = math.Pi / 180
	rad2deg = 180 / math.Pi
)

// SkyCoord stores sky coordinates.
type SkyCoord struct {
	RA  float64 // Right ascension (or galactic/ecliptic longitude)
	Dec float64 // Declination (or latitude)
}

// Pixel represents a single pixel in form of slice int{x, y}.
type Pixel [2]int // x, y

// WCS represents a World Coordinate System transformation defined by a FITS header.
// It maps 2D detector pixel coordinates to celestial sky coordinates (such as ICRS
// Right Ascension and Declination) and vice versa, handling linear scaling, rotation,
// skew, and spherical projections.
//
// Pixel coordinates use 0-based indexing, where (0, 0) corresponds to the center of
// the bottom-left pixel (equivalent to (1, 1) in standard FITS convention).
//
// A WCS instance is safe for concurrent read-only transformations across multiple goroutines.
type WCS struct {
	crpix [2]float64    // reference pixel in FITS 1-based convention (CRPIX1, CRPIX2)
	crval [2]float64    // reference world coordinates in degrees (CRVAL1, CRVAL2)
	cd    [2][2]float64 // linear transformation in degrees per pixel: CD, or CDELT * PC
	cdInv [2][2]float64 // inverse of cd, used for world-to-pixel conversion
	naxis [2]int        // image size in pixels: NAXIS1 (width), NAXIS2 (height)
}

// NewFromHeaders creates a WCS from the header of an image HDU (primary or IMAGE extension).
//
// The header must contain NAXIS1 and NAXIS2 (image size) and CTYPE1 and CTYPE2 with
// the gnomonic projection code "TAN" (e.g. "RA---TAN", "DEC--TAN"). CRPIX1/2 and CRVAL1/2
// default to 0 when absent, as defined by the FITS WCS standard.
//
// The linear transformation is taken from the CDi_j keywords when any of them is present.
// Otherwise it is built as CDELTi * PCi_j, where PCi_j defaults to the identity matrix
// and CDELTi defaults to 1. Units must be degrees (CUNITi absent or "deg").
//
// An error is returned for missing or invalid keywords, unsupported projections or units,
// and a singular (non-invertible) transformation matrix.
func NewFromHeaders(h fits.Headers) (*WCS, error) {
	var w WCS
	var err error

	if w.naxis[0], err = headerInt(h, "NAXIS1"); err != nil {
		return nil, err
	}
	if w.naxis[1], err = headerInt(h, "NAXIS2"); err != nil {
		return nil, err
	}
	if w.naxis[0] <= 0 || w.naxis[1] <= 0 {
		return nil, fmt.Errorf("invalid image size %dx%d", w.naxis[0], w.naxis[1])
	}

	for i := range 2 {
		if w.crpix[i], err = headerFloat(h, fmt.Sprintf("CRPIX%d", i+1), 0); err != nil {
			return nil, err
		}
		if w.crval[i], err = headerFloat(h, fmt.Sprintf("CRVAL%d", i+1), 0); err != nil {
			return nil, err
		}

		ctype, err := headerString(h, fmt.Sprintf("CTYPE%d", i+1))
		if err != nil {
			return nil, err
		}

		// CTYPE is 8 characters: "XXXX-TAN", where the projection code occupies positions 5-7.
		if len(ctype) != 8 || ctype[4] != '-' || ctype[5:] != "TAN" {
			return nil, fmt.Errorf("unsupported projection %q: only TAN is supported", ctype)
		}

		unit := strings.TrimSpace(h[fmt.Sprintf("CUNIT%d", i+1)])
		if unit != "" && !strings.EqualFold(unit, "deg") {
			return nil, fmt.Errorf("unsupported unit %q: only deg is supported", unit)
		}
	}

	cd, err := parseMatrix(h)
	if err != nil {
		return nil, err
	}

	det := cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0]
	if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
		return nil, fmt.Errorf("singular transformation matrix: determinant is %g", det)
	}

	w.cd = cd
	w.cdInv = [2][2]float64{
		{cd[1][1] / det, -cd[0][1] / det},
		{-cd[1][0] / det, cd[0][0] / det},
	}

	return &w, nil
}

// parseMatrix returns the linear transformation matrix in degrees per pixel.
// If any CDi_j keyword is present, the CD matrix is used and missing elements are 0.
// Otherwise the matrix is CDELTi * PCi_j with identity and unit defaults.
func parseMatrix(h fits.Headers) ([2][2]float64, error) {
	var m [2][2]float64

	hasCD := false
	for i := range 2 {
		for j := range 2 {
			if _, ok := h[fmt.Sprintf("CD%d_%d", i+1, j+1)]; ok {
				hasCD = true
			}
		}
	}

	if hasCD {
		for i := range 2 {
			for j := range 2 {
				v, err := headerFloat(h, fmt.Sprintf("CD%d_%d", i+1, j+1), 0)
				if err != nil {
					return m, err
				}

				m[i][j] = v
			}
		}

		return m, nil
	}

	var cdelt [2]float64
	for i := range 2 {
		v, err := headerFloat(h, fmt.Sprintf("CDELT%d", i+1), 1)
		if err != nil {
			return m, err
		}

		cdelt[i] = v
	}

	for i := range 2 {
		for j := range 2 {
			def := 0.0
			if i == j {
				def = 1
			}

			pc, err := headerFloat(h, fmt.Sprintf("PC%d_%d", i+1, j+1), def)
			if err != nil {
				return m, err
			}

			m[i][j] = cdelt[i] * pc
		}
	}

	return m, nil
}

// headerFloat returns the floating point value of the key, or def if the key is absent.
// FITS allows the 'D' exponent marker (e.g. 1.0D-03), which is converted to 'E'.
// NaN and infinite values are rejected.
func headerFloat(h fits.Headers, key string, def float64) (float64, error) {
	raw, ok := h[key]
	if !ok {
		return def, nil
	}

	s := strings.NewReplacer("D", "E", "d", "E").Replace(strings.TrimSpace(raw))

	v, err := strconv.ParseFloat(s, 64)
	if err != nil || math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, fmt.Errorf("invalid %s: %q", key, raw)
	}

	return v, nil
}

// headerInt returns the integer value of a required key.
func headerInt(h fits.Headers, key string) (int, error) {
	raw, ok := h[key]
	if !ok {
		return 0, fmt.Errorf("missing %s", key)
	}

	v, err := strconv.Atoi(strings.TrimSpace(raw))
	if err != nil {
		return 0, fmt.Errorf("invalid %s: %q", key, raw)
	}

	return v, nil
}

// headerString returns the trimmed value of a required key.
func headerString(h fits.Headers, key string) (string, error) {
	raw, ok := h[key]
	if !ok {
		return "", fmt.Errorf("missing %s", key)
	}

	return strings.TrimSpace(raw), nil
}

// PixelToWorld converts a single point from an image into celestial coordinates.
//
// The pixel uses 0-based indexing (see WCS). The returned coordinates are in degrees.
// The TAN projection is defined for every pixel, so the error is always nil.
func (w *WCS) PixelToWorld(p Pixel) (*SkyCoord, error) {
	c := w.pixelToWorld(float64(p[0]), float64(p[1]))
	return &c, nil
}

// PixelsToWorld converts a slice of points from an image into slice of celestial coordinates.
//
// The result has the same length and order as the input. Pixels use 0-based indexing,
// and the returned coordinates are in degrees. The error is always nil for the TAN projection.
func (w *WCS) PixelsToWorld(pixels []Pixel) ([]SkyCoord, error) {
	coords := make([]SkyCoord, len(pixels))
	for i, p := range pixels {
		coords[i] = w.pixelToWorld(float64(p[0]), float64(p[1]))
	}

	return coords, nil
}

// WorldToPixel converts celestial coordinates into a pixel position on an image.
//
// The returned pixel is the one whose area contains the exact (fractional) position,
// i.e. the nearest pixel center. Coordinates are in degrees; declination must be within [-90, 90].
// An error is returned if the coordinate lies at or beyond 90 degrees from the reference point
// (CRVAL), where the TAN projection is undefined.
func (w *WCS) WorldToPixel(coord SkyCoord) (Pixel, error) {
	x, y, err := w.worldToPixel(coord)
	if err != nil {
		return Pixel{}, err
	}

	return Pixel{int(math.Floor(x + 0.5)), int(math.Floor(y + 0.5))}, nil
}

// WorldToPixels converts a slice of celestial coordinates into a slice of pixels position on an image.
//
// The result has the same length and order as the input. The conversion stops at the first
// coordinate that cannot be projected, and the error includes its index.
func (w *WCS) WorldToPixels(coords []SkyCoord) ([]Pixel, error) {
	pixels := make([]Pixel, len(coords))
	for i, c := range coords {
		p, err := w.WorldToPixel(c)
		if err != nil {
			return nil, fmt.Errorf("coordinate #%d: %w", i, err)
		}

		pixels[i] = p
	}

	return pixels, nil
}

// PixelScale returns the image scale in arcseconds per pixel (arcsec/pixel) calculated from the CD or CDELT matrix.
//
// The value is the geometric mean of the two axis scales, sqrt(|det(CD)|), which is exact for
// square pixels and a common approximation otherwise. The error is always nil.
func (w *WCS) PixelScale() (float64, error) {
	det := w.cd[0][0]*w.cd[1][1] - w.cd[0][1]*w.cd[1][0]
	return math.Sqrt(math.Abs(det)) * 3600, nil
}

// Footprint returns the sky coordinates for the 4 corners of the image.
//
// The corners are the outer edges of the first and last pixels, in order:
// bottom-left, bottom-right, top-right, top-left. Coordinates are in degrees.
func (w *WCS) Footprint() []SkyCoord {
	x0, y0 := -0.5, -0.5
	x1 := float64(w.naxis[0]) - 0.5
	y1 := float64(w.naxis[1]) - 0.5

	corners := [4][2]float64{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}
	coords := make([]SkyCoord, len(corners))
	for i, c := range corners {
		coords[i] = w.pixelToWorld(c[0], c[1])
	}

	return coords
}

// Contains returns true if the given celestial object lies within the boundaries of the matrix.
//
// The boundaries are the outer edges of the image (pixels 0..NAXIS1-1 and 0..NAXIS2-1 in
// 0-based indexing). It returns false for coordinates that cannot be projected.
func (w *WCS) Contains(coord SkyCoord) bool {
	x, y, err := w.worldToPixel(coord)
	if err != nil {
		return false
	}

	px := math.Floor(x + 0.5)
	py := math.Floor(y + 0.5)

	return px >= 0 && px < float64(w.naxis[0]) && py >= 0 && py < float64(w.naxis[1])
}

// Slice returns new WCS object with an updated reference point (CRPIX) after cropping the image.
//
// The crop region is half-open, as in Go slicing: columns [xMin, xMax) and rows [yMin, yMax),
// 0-based. The region must lie within the image and be non-empty. The original WCS is not modified,
// and the new WCS describes the cropped image, whose pixel (0, 0) is the original pixel (xMin, yMin).
func (w *WCS) Slice(xMin, yMin, xMax, yMax int) (*WCS, error) {
	if xMin < 0 || yMin < 0 || xMin >= xMax || yMin >= yMax || xMax > w.naxis[0] || yMax > w.naxis[1] {
		return nil, fmt.Errorf(
			"invalid slice [%d:%d, %d:%d] for image of size %dx%d",
			xMin, xMax, yMin, yMax, w.naxis[0], w.naxis[1],
		)
	}

	s := *w
	s.crpix[0] -= float64(xMin)
	s.crpix[1] -= float64(yMin)
	s.naxis = [2]int{xMax - xMin, yMax - yMin}

	return &s, nil
}

// pixelToWorld converts 0-based (possibly fractional) pixel coordinates into sky coordinates in degrees.
// The conversion is linear (CD matrix) followed by the inverse gnomonic (TAN) projection.
// RA is normalized to [0, 360).
func (w *WCS) pixelToWorld(x, y float64) SkyCoord {
	// 0-based -> 1-based, relative to the reference pixel.
	dx := x + 1 - w.crpix[0]
	dy := y + 1 - w.crpix[1]

	// Intermediate world coordinates in rad.
	xi := (w.cd[0][0]*dx + w.cd[0][1]*dy) * deg2rad
	eta := (w.cd[1][0]*dx + w.cd[1][1]*dy) * deg2rad

	// TAN deprojection around the reference point (CRVAL).
	alpha0 := w.crval[0] * deg2rad
	delta0 := w.crval[1] * deg2rad
	sin0, cos0 := math.Sincos(delta0)

	d := cos0 - eta*sin0
	alpha := alpha0 + math.Atan2(xi, d)
	delta := math.Atan2(eta*cos0+sin0, math.Hypot(xi, d))

	ra := math.Mod(alpha*rad2deg, 360)
	if ra < 0 {
		ra += 360
	}

	return SkyCoord{RA: ra, Dec: delta * rad2deg}
}

// worldToPixel converts sky coordinates in degrees into 0-based (possibly fractional) pixel coordinates.
// It applies the gnomonic (TAN) projection around the reference point (CRVAL), followed by the
// inverse linear transformation (CD matrix). It returns an error if the declination is out of range
// or the point is at or beyond 90 degrees from the reference point.
func (w *WCS) worldToPixel(c SkyCoord) (float64, float64, error) {
	if math.IsNaN(c.RA) || math.IsInf(c.RA, 0) || !(math.Abs(c.Dec) <= 90) {
		return 0, 0, fmt.Errorf("invalid coordinate RA=%g Dec=%g: RA must be finite and declination must be within [-90, 90]", c.RA, c.Dec)
	}

	alpha0 := w.crval[0] * deg2rad
	delta0 := w.crval[1] * deg2rad
	sin0, cos0 := math.Sincos(delta0)

	sinD, cosD := math.Sincos(c.Dec * deg2rad)
	sinA, cosA := math.Sincos(c.RA*deg2rad - alpha0)

	cosC := sinD*sin0 + cosD*cos0*cosA
	if cosC <= 0 {
		return 0, 0, fmt.Errorf("coordinate RA=%g Dec=%g is at or beyond 90 degrees from the reference point", c.RA, c.Dec)
	}

	// Standard world coordinates in degrees.
	xi := (cosD * sinA / cosC) * rad2deg
	eta := ((sinD*cos0 - cosD*sin0*cosA) / cosC) * rad2deg

	// Inverse linear transformation.
	dx := w.cdInv[0][0]*xi + w.cdInv[0][1]*eta
	dy := w.cdInv[1][0]*xi + w.cdInv[1][1]*eta

	// 1-based -> 0-based.
	return dx + w.crpix[0] - 1, dy + w.crpix[1] - 1, nil
}
