package fits

// ExtensionType is the XTENSION keyword value identifying the kind of a FITS
// extension HDU.
type ExtensionType string

const (
	// ExtImage denotes an image extension (XTENSION = 'IMAGE').
	ExtImage ExtensionType = "IMAGE"
	// ExtTable denotes an ASCII table extension (XTENSION = 'TABLE').
	ExtTable ExtensionType = "TABLE"
	// ExtBinTable denotes a binary table extension (XTENSION = 'BINTABLE').
	ExtBinTable ExtensionType = "BINTABLE"
)
