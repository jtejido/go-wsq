# go-wsq
Golang port of NBIS' WSQ Image Encoder and Decoder.

## Design
It follows Go's standard decoder and encoder (jpeg, png, gif, bmp) for your nifty interfaces.

### Encode
```
Encode(w io.Writer, img image.Image, o *Options) error
```

The options are the following:
```
type Options struct {
	Bitrate  float32
	Comments []string
	Metadata map[string]string
}
```
### Decode
```
Decode(r io.Reader) (image.Image, error)
```
The decode API will possibly change to return the metadata and comments added, at the moment it is just kept in memory and discarded.

### Usage
See the sample folder for the simple **decode-encode-decode-convert** test.

### Objective
This has been created for the Golang port of SourceAFIS.