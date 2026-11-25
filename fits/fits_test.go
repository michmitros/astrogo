package fits

import "testing"

func TestHDUIsPrimary(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    bool
	}{
		"SIMPLE only": {
			headers: Headers{"SIMPLE": "T"},
			want:    true,
		},
		"SIMPLE and XTENSION": {
			headers: Headers{"SIMPLE": "T", "XTENSION": "IMAGE"},
			want:    false,
		},
		"XTENSION only": {
			headers: Headers{"XTENSION": "IMAGE"},
			want:    false,
		},
		"no keywords": {
			headers: Headers{},
			want:    false,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.IsPrimary(); got != tt.want {
				t.Errorf("IsPrimary() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestHDUBitpix(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    int
	}{
		"8 bits": {
			headers: Headers{"BITPIX": "8"},
			want:    8,
		},
		"negative 32 bits": {
			headers: Headers{"BITPIX": "-32"},
			want:    -32,
		},
		"missing keyword": {
			headers: Headers{},
			want:    0,
		},
		"invalid value": {
			headers: Headers{"BITPIX": "abc"},
			want:    0,
		},
		"value out of int32 range": {
			headers: Headers{"BITPIX": "9999999999"},
			want:    0,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.Bitpix(); got != tt.want {
				t.Errorf("Bitpix() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestHDUNaxis(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    int
	}{
		"two axes": {
			headers: Headers{"NAXIS": "2"},
			want:    2,
		},
		"zero axes": {
			headers: Headers{"NAXIS": "0"},
			want:    0,
		},
		"missing keyword": {
			headers: Headers{},
			want:    0,
		},
		"invalid value": {
			headers: Headers{"NAXIS": "two"},
			want:    0,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.Naxis(); got != tt.want {
				t.Errorf("Naxis() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestHDUNaxisN(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		n       uint
		want    int
	}{
		"first axis": {
			headers: Headers{"NAXIS1": "640"},
			n:       1,
			want:    640,
		},
		"second axis": {
			headers: Headers{"NAXIS1": "640", "NAXIS2": "480"},
			n:       2,
			want:    480,
		},
		"missing axis": {
			headers: Headers{"NAXIS1": "640"},
			n:       2,
			want:    0,
		},
		"invalid axis value": {
			headers: Headers{"NAXIS1": "x"},
			n:       1,
			want:    0,
		},
		"axis zero is not stored": {
			headers: Headers{"NAXIS1": "640"},
			n:       0,
			want:    0,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.NaxisN(tt.n); got != tt.want {
				t.Errorf("NaxisN(%d) = %d, want %d", tt.n, got, tt.want)
			}
		})
	}
}

func TestHDUXtension(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    string
	}{
		"image extension": {
			headers: Headers{"XTENSION": "IMAGE"},
			want:    "IMAGE",
		},
		"binary table": {
			headers: Headers{"XTENSION": "BINTABLE"},
			want:    "BINTABLE",
		},
		"keyword absent": {
			headers: Headers{"SIMPLE": "T"},
			want:    "",
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.Extension(); string(got) != tt.want {
				t.Errorf("Xtension() = %q, want %q", got, tt.want)
			}
		})
	}
}

func TestHDUExtName(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    string
	}{
		"extension name present": {
			headers: Headers{"EXTNAME": "SCI"},
			want:    "SCI",
		},
		"keyword absent": {
			headers: Headers{},
			want:    "",
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.ExtName(); got != tt.want {
				t.Errorf("ExtName() = %q, want %q", got, tt.want)
			}
		})
	}
}

func TestHDUIsImage(t *testing.T) {
	tests := map[string]struct {
		headers Headers
		want    bool
	}{
		"exact IMAGE": {
			headers: Headers{"XTENSION": "IMAGE"},
			want:    true,
		},
		"BINTABLE is not an image": {
			headers: Headers{"XTENSION": "BINTABLE"},
			want:    false,
		},
		"lowercase image is not an image": {
			headers: Headers{"XTENSION": "image"},
			want:    false,
		},
		"primary HDU is not an image": {
			headers: Headers{"SIMPLE": "T"},
			want:    false,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			hdu := HDU{Headers: tt.headers}
			if got := hdu.IsImage(); got != tt.want {
				t.Errorf("IsImage() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestParseIntHeader(t *testing.T) {
	tests := map[string]struct {
		raw  string
		want int
	}{
		"positive number": {
			raw:  "123",
			want: 123,
		},
		"negative number": {
			raw:  "-64",
			want: -64,
		},
		"empty string": {
			raw:  "",
			want: 0,
		},
		"not a number": {
			raw:  "1.5",
			want: 0,
		},
		"exceeds int32 range": {
			raw:  "2147483648",
			want: 0,
		},
	}

	for name, tt := range tests {
		t.Run(name, func(t *testing.T) {
			if got := parseIntHeader(tt.raw); got != tt.want {
				t.Errorf("parseIntHeader(%q) = %d, want %d", tt.raw, got, tt.want)
			}
		})
	}
}
