[![PkgGoDev](https://pkg.go.dev/badge/github.com/michmitros/astrogo/fits)](https://pkg.go.dev/github.com/michmitros/astrogo/fits)

# fits

The `fits` package reads files in the [FITS](https://fits.gsfc.nasa.gov/) (Flexible Image Transport System) format, the standard data format in astronomy. Currently the package supports **reading** files only.

```bash
go get github.com/michmitros/astrogo/fits
```

## FITS file structure

A FITS file consists of Header Data Units (HDUs). Each HDU contains a header (a set of 80-byte cards) and optional data. The first HDU is always the primary HDU (`SIMPLE`), followed by optional extensions (`XTENSION`).

The package uses constants defined by the FITS standard:

- `BLOCK_SIZE` = 2880 bytes – the logical block length; headers and data are padded to a multiple of this size,
- `HDU_LINE_LENGTH` = 80 bytes – the length of a single header card.

## Types

### `FitsFile`

```go
type FitsFile []HDU
```

The result of decoding a whole file: a slice of HDUs in the order they appear in the file. Element `[0]` is always the primary HDU.

### `HDU`

```go
type HDU struct {
    Headers Headers
    Data    []byte
}
```

- `Headers` – the parsed header cards,
- `Data` – raw data bytes, without padding. For `NAXIS = 0` the field is empty. Multi-byte data types (`BITPIX` 16, 32, -32, -64) are already converted to **Little Endian** byte order.

Methods of `HDU`:

| Method | Description |
|---|---|
| `IsPrimary() bool` | Returns `true` if the header contains `SIMPLE` and does not contain `XTENSION`. |
| `Bitpix() int` | The `BITPIX` value (e.g. 8, 16, 32, -32, -64). Returns `0` if the keyword is missing or is not a valid integer. |
| `Naxis() int` | The `NAXIS` value, i.e. the number of data axes. Returns `0` if missing or invalid. |
| `NaxisN(n uint) int` | The length of the n-th data axis (`NAXISn`), numbered from 1. Returns `0` if missing or invalid. |
| `Xtension() string` | The `XTENSION` value (e.g. `"IMAGE"`, `"BINTABLE"`). Empty string if the keyword is absent. |
| `ExtName() string` | The `EXTNAME` value. Empty string if the keyword is absent. |
| `IsImage() bool` | Returns `true` if `XTENSION` is exactly `"IMAGE"`. Primary HDUs are not reported as images. |

### `Headers`

```go
type Headers map[string]string
```

A map of header cards:

- keys are card names without spaces (e.g. `"NAXIS1"`),
- values are raw strings: quotes are removed, doubled apostrophes (`''`) are replaced by a single apostrophe, and the comment after `/` is dropped,
- `COMMENT`, `HISTORY`, blank and `END` cards are not stored; cards without `=` in column 9 are skipped as well,
- if a key appears more than once, the last value is kept.

Because cards without `=` are skipped, long string values (`CONTINUE`) and the `HIERARCH` convention are not supported.

## Reading a file

### `Decode`

```go
func Decode(r io.Reader) (FitsFile, error)
```

Reads the whole file from `r` and returns all HDUs in order. Behavior:

- the first HDU must be primary (`SIMPLE`), otherwise an error is returned,
- if `r` contains no HDU (an empty file or an empty block at the beginning), `ErrEmptyFile` is returned,
- errors from a specific HDU are wrapped with its 0-based index.

### `Decoder`

```go
func NewDecoder(r io.Reader) *Decoder
func (d *Decoder) HDU() (HDU, error)
```

Reads HDUs from a stream one at a time. Unlike `Decode`, it does not buffer the whole file, so it suits streams with many HDUs.

- If the stream is empty, `ErrEmptyFile` is returned.
- The check that the first HDU is primary is **not performed**.
- After the last HDU, the next call returns an error describing the end of the stream (`failed to decode HDU #N: EOF`). Because the error is formatted with `%s` instead of `%w`, it **cannot be matched with `errors.Is(err, io.EOF)`**. A reliable way to stop reading is to know the number of HDUs in advance or to stop at the first error.

### Errors

- `ErrEmptyFile` – no HDU was found,
- validation errors (e.g. missing `NAXIS`, invalid `BITPIX`, `NAXISn <= 0`, data size exceeding the `2^60` byte limit, invalid first header key) – returned with a description of the problem,
- premature end of data – `io.ErrUnexpectedEOF` is returned.

Supported `BITPIX` values: `8`, `16`, `32`, `-32`, `-64`. The maximum `NAXIS` value is 999.

## Image conversion

### `ToImage`

```go
func ToImage(hdu HDU, zscaleContrast float64, zscaleSamples int) (*image.Gray16, error)
```

Converts a two-dimensional image HDU (`NAXIS = 2`) into a 16-bit grayscale image (`*image.Gray16` from the standard `image` package):

1. Pixel values are computed as `BZERO + BSCALE * value` (defaults `BZERO = 0`, `BSCALE = 1`).
2. The display range is determined with the **ZScale** algorithm:
   - `zscaleContrast` – the contrast parameter; larger values give a narrower range,
   - `zscaleSamples` – the approximate maximum number of sampled pixels; a value `<= 0` means all pixels are used.
3. Values within the range are mapped non-linearly to `0..65535` using `2/π · asin(t)`. Values below the range become 0, values above it become the maximum.
4. `NaN` and infinite pixels are rendered as 0 (black). For an empty range (e.g. a constant image) all pixels are 0.
5. The image is flipped vertically: FITS stores the first row at the bottom, while the result has it at the top.

An error is returned when `NAXIS != 2`, when `NAXIS1` or `NAXIS2` are not positive, or when there is not enough data for `NAXIS1 * NAXIS2` pixels of `BITPIX` size. Unsupported `BITPIX` values are treated as zeros.

## Examples

Reading a file, selecting the primary HDU and saving the image as PNG:

```go
package main

import (
    "image/png"
    "log"
    "os"

    "github.com/michmitros/astrogo/fits"
)

func main() {
    f, err := os.Open("image.fits")
    if err != nil {
        log.Fatal(err)
    }
    defer f.Close()

    file, err := fits.Decode(f)
    if err != nil {
        log.Fatal(err)
    }

    img, err := fits.ToImage(file[0], 0.25, 1000)
    if err != nil {
        log.Fatal(err)
    }

    out, err := os.Create("image.png")
    if err != nil {
        log.Fatal(err)
    }
    defer out.Close()

    if err := png.Encode(out, img); err != nil {
        log.Fatal(err)
    }
}
```

Streaming two HDUs (the number of HDUs is known in advance):

```go
d := fits.NewDecoder(f)
for i := 0; i < 2; i++ {
    hdu, err := d.HDU()
    if err != nil {
        log.Fatal(err)
    }
    log.Printf("HDU: XTENSION=%q EXTNAME=%q NAXIS=%d", hdu.Xtension(), hdu.ExtName(), hdu.Naxis())
}
```

## Limitations

- Read-only; FITS files cannot be written.
- Tables (`BINTABLE`, `TABLE`) are parsed with `ToTable`; `TSCALn`/`TZEROn` scaling is not applied to column data (the values are exposed in `ColumnMeta`).
- `ToImage` supports only two-dimensional images.
- `CONTINUE` and `HIERARCH` cards are not supported.
- Compression and the group data convention (`GROUPS`) are not supported.
