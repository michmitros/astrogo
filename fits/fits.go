package fits

import (
	"fmt"
	"strconv"
)

// FitsFile is a slice of HDUs from single file. The first element is always
// the primary HDU, followed by any extensions in the order they appear in the file.
type FitsFile []HDU

// Headers is a map of HDU headers. Keys are card names with surrounding spaces
// removed (e.g. "NAXIS1"). Values are the raw card values as strings: quoted
// strings are unquoted (with ” unescaped to '), and the "/" comment is dropped.
// COMMENT, HISTORY, blank and END cards are not stored, and neither are cards
// without "=" in column 9. If a key appears more than once, the last card wins.
type Headers map[string]string

// HDU (Header Data Unit) is the highest level component of the FITS file structure,
// consisting of a header and usually a data. Data stored in HDU (if present) are
// already converted to Little Endian.
type HDU struct {
	// Headers contains the parsed header cards of the HDU.
	Headers Headers
	// Data contains the raw data bytes of the HDU, padded data excluded.
	// It is empty when NAXIS is 0.
	Data []byte
}

// IsPrimary reports whether the HDU is a primary HDU, i.e. its header contains
// the SIMPLE keyword and does not contain XTENSION.
func (u HDU) IsPrimary() bool {
	return isPrimary(u.Headers)
}

// Bitpix returns the BITPIX value of the HDU (bits per pixel, e.g. 8, 16, 32,
// -32 or -64). It returns 0 if BITPIX is missing or not a valid integer.
func (u HDU) Bitpix() int {
	return bitpix(u.Headers)
}

// Naxis returns the NAXIS value of the HDU, i.e. the number of data axes.
// It returns 0 if NAXIS is missing or not a valid integer.
func (u HDU) Naxis() int {
	return naxis(u.Headers)
}

// NaxisN returns the length of the n-th data axis (NAXISn) of the HDU.
// Axes are numbered from 1, so NaxisN(1) is the NAXIS1 value.
// It returns 0 if NAXISn is missing or not a valid integer.
func (u HDU) NaxisN(n uint) int {
	return naxisN(u.Headers, n)
}

// Xtension returns the XTENSION value of the HDU (e.g. "IMAGE", "BINTABLE").
// It returns an empty string if the keyword is absent, as in a primary HDU.
func (u HDU) Extension() ExtensionType {
	return ExtensionType(u.Headers["XTENSION"])
}

// ExtName returns the EXTNAME value of the HDU. It returns an empty string
// if the keyword is absent.
func (u HDU) ExtName() string {
	return u.Headers["EXTNAME"]
}

// IsImage reports whether the HDU is an image extension, i.e. its XTENSION
// value is exactly "IMAGE". Primary HDUs are not reported as images.
func (u HDU) IsImage() bool {
	return u.Headers["XTENSION"] == "IMAGE"
}

func isPrimary(h Headers) bool {
	_, hasSimple := h["SIMPLE"]
	_, hasXtension := h["XTENSION"]
	return hasSimple && !hasXtension
}

func bitpix(h Headers) int {
	return parseIntHeader(h["BITPIX"])
}

func naxis(h Headers) int {
	return parseIntHeader(h["NAXIS"])
}

func naxisN(h Headers, n uint) int {
	return parseIntHeader(h[fmt.Sprintf("NAXIS%d", n)])
}

func parseIntHeader(h string) int {
	if val, err := strconv.ParseInt(h, 10, 32); err == nil {
		return int(val)
	}
	return 0
}
