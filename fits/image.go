package fits

import (
	"encoding/binary"
	"fmt"
	"image"
	"image/color"
	"math"
	"slices"
	"strconv"
)

// ToImage converts a 2D image HDU into a 16-bit grayscale image.
//
// Raw pixel values are scaled with BZERO and BSCALE header keywords (defaults 0 and 1).
// The display range is computed with the ZScale algorithm (see zscale): zscaleContrast
// is the contrast parameter (larger values give a narrower range) and zscaleSamples
// is the approximate maximum number of pixels sampled for the estimate; values <= 0
// mean that all pixels are used. Values inside the range are mapped non-linearly to
// 0..65535 using 2/pi*asin(t). Values below the range become 0 and values above it
// become the maximum. NaN and infinite pixels are rendered as 0 (black). If the range
// is empty (for example, a constant image), all pixels are rendered as 0.
//
// The image is flipped vertically: FITS stores the first row at the bottom, while
// the returned image has it at the top.
//
// An error is returned if NAXIS is not 2, if NAXIS1 or NAXIS2 are not positive,
// or if the data is shorter than required by NAXIS1*NAXIS2 pixels of BITPIX size.
// Unsupported BITPIX values are read as zero.
func ToImage(hdu HDU, zscaleContrast float64, zscaleSamples int) (*image.Gray16, error) {
	if ext := hdu.Extension(); ext != ExtImage {
		return nil, fmt.Errorf("hdu is not an image, got XTENSION=%s", ext)
	}

	if naxis := hdu.Naxis(); naxis != 2 {
		return nil, fmt.Errorf("expected 2D image, got NAXIS=%d", naxis)
	}

	w, h := hdu.NaxisN(1), hdu.NaxisN(2)
	if w <= 0 || h <= 0 {
		return nil, fmt.Errorf("invalid image size %dx%d", w, h)
	}

	bitpix := hdu.Bitpix()
	bytesPerPixel := bitpix
	if bytesPerPixel < 0 {
		bytesPerPixel = -bytesPerPixel
	}
	bytesPerPixel /= 8

	if len(hdu.Data) < w*h*bytesPerPixel {
		return nil, fmt.Errorf("data too short: have %d bytes, need %d", len(hdu.Data), w*h*bytesPerPixel)
	}

	bzero := floatHeader(hdu.Headers, "BZERO", 0)
	bscale := floatHeader(hdu.Headers, "BSCALE", 1)

	values := make([]float64, w*h)
	for i := range values {
		values[i] = bzero + bscale*readValue(hdu.Data, bitpix, i)
	}

	z1, z2 := zscale(values, zscaleContrast, zscaleSamples)
	span := z2 - z1

	img := image.NewGray16(image.Rect(0, 0, w, h))
	for i, v := range values {
		var g uint16
		if span > 0 && !math.IsNaN(v) && !math.IsInf(v, 0) {
			t := min(max((v-z1)/span, 0), 1)
			g = uint16(65535 * 2 / math.Pi * math.Asin(t))
		}

		x, y := i%w, i/w
		img.Set(x, h-1-y, color.Gray16{Y: g})
	}

	return img, nil
}

func zscale(values []float64, contrast float64, numOfSamples int) (float64, float64) {
	stride := 1
	if numOfSamples > 0 && len(values) > numOfSamples {
		stride = len(values) / numOfSamples
	}

	samples := make([]float64, 0, min(len(values), max(numOfSamples, 0)))
	for i := 0; i < len(values); i += stride {
		v := values[i]
		if math.IsNaN(v) || math.IsInf(v, 0) {
			continue
		}
		samples = append(samples, v)
	}

	if len(samples) == 0 {
		return 0, 0
	}

	slices.Sort(samples)

	const (
		minNPixels    = 5
		krej          = 2.5
		maxIterations = 5
		maxReject     = 0.5
	)

	npix := len(samples)
	lo, hi := samples[0], samples[npix-1]
	center := npix / 2
	median := samples[center]
	minPix := max(minNPixels, int(float64(npix)*maxReject))

	good := make([]bool, npix)
	for i := range good {
		good[i] = true
	}

	ngood := npix
	lastNgood := npix + 1
	var slope, intercept float64

	for iter := 0; iter < maxIterations; iter++ {
		if ngood >= lastNgood || ngood < minPix {
			break
		}

		var sx, sy, sxx, sxy, n float64
		for i, y := range samples {
			if !good[i] {
				continue
			}
			x := float64(i)
			sx += x
			sy += y
			sxx += x * x
			sxy += x * y
			n++
		}

		denom := n*sxx - sx*sx
		if denom == 0 {
			break
		}
		slope = (n*sxy - sx*sy) / denom
		intercept = (sy - slope*sx) / n

		var sr, srr, nr float64
		for i, y := range samples {
			if !good[i] {
				continue
			}
			r := y - (intercept + slope*float64(i))
			sr += r
			srr += r * r
			nr++
		}
		mean := sr / nr
		sigma := math.Sqrt(max(srr/nr-mean*mean, 0))
		threshold := krej * sigma

		lastNgood = ngood
		ngood = 0
		for i, y := range samples {
			r := y - (intercept + slope*float64(i))
			good[i] = math.Abs(r) <= threshold
			if good[i] {
				ngood++
			}
		}
	}

	if ngood < minPix {
		return lo, hi
	}

	slope /= contrast
	z1 := max(lo, median-float64(center)*slope)
	z2 := min(hi, median+float64(npix-1-center)*slope)

	return z1, z2
}

func readValue(data []byte, bitpix, i int) float64 {
	switch bitpix {
	case 8:
		return float64(data[i])
	case 16:
		return float64(int16(binary.LittleEndian.Uint16(data[2*i:])))
	case 32:
		return float64(int32(binary.LittleEndian.Uint32(data[4*i:])))
	case -32:
		return float64(math.Float32frombits(binary.LittleEndian.Uint32(data[4*i:])))
	case -64:
		return math.Float64frombits(binary.LittleEndian.Uint64(data[8*i:]))
	}
	return 0
}

func floatHeader(h Headers, key string, def float64) float64 {
	raw, ok := h[key]
	if !ok {
		return def
	}

	v, err := strconv.ParseFloat(raw, 64)
	if err != nil {
		return def
	}

	return v
}
