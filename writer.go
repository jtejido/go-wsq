package wsq

import (
	"bufio"
	"fmt"
	"image"
	"io"
	"log"
	"math"
	"net/url"
	"strings"
)

const DefaultBitrate float32 = .75

type Options struct {
	Bitrate  float32
	Comments []string
	Metadata map[string]string
}

func Encode(w io.Writer, img image.Image, o *Options) error {
	if img == nil {
		return fmt.Errorf("image cannot be nil.")
	}
	var m *image.Gray
	switch img.(type) {
	case *image.Gray:
		m = img.(*image.Gray)
	default:
		return fmt.Errorf("only supports image.Gray at the moment.")
	}
	return encode(w, m, o)
}

func encode(w io.Writer, m *image.Gray, o *Options) (err error) {
	b := m.Bounds()
	var e encoder
	if ww, ok := w.(writer); ok {
		e.w = ww
	} else {
		e.w = bufio.NewWriter(w)
	}
	var fdata []float32 /* floating point pixel image  */
	var qdata []int     /* quantized image pointer     */
	/* quantized block sizes */
	qsize := &reference[int]{0}
	qsize1 := &reference[int]{0}
	qsize2 := &reference[int]{0}
	qsize3 := &reference[int]{0}
	hufftable := []huffCode{}
	huffbits := &reference[[]int]{}
	huffvalues := &reference[[]int]{} /* huffman code parameters */

	/* Convert image pixels to floating point. */
	m_shift := &reference[float32]{0}
	r_scale := &reference[float32]{0}
	fdata = convertImageToFloat(m, m_shift, r_scale)
	e.init()

	/* Build WSQ decomposition trees */
	e.buildWSQTrees(b.Dx(), b.Dy())

	/* WSQ decompose the image */
	err = wsqDecompose(&e, fdata, b.Dx(), b.Dy(), e.tableDTT.hifilt, max_hifilt, e.tableDTT.lofilt, max_lofilt)
	if err != nil {
		return
	}
	/* Set compression ratio and 'q' to zero. */
	e.quant_vals.cr = 0
	e.quant_vals.q = 0.0

	/* Assign specified r-bitrate into quantization structure. */
	e.quant_vals.r = DefaultBitrate
	if o != nil {
		e.quant_vals.r = o.Bitrate
		if e.quant_vals.r <= 0 {
			e.quant_vals.r = DefaultBitrate
		}
	}

	/* Compute subband variances. */
	variance(&e, fdata, b.Dx(), b.Dy())

	/* Quantize the floating point pixmap. */

	qdata = quantize(&e, qsize, fdata, b.Dx(), b.Dy())

	/* Compute quantized WSQ subband block sizes */
	quantBlockSizes(&e, qsize1, qsize2, qsize3)

	if qsize.value != qsize1.value+qsize2.value+qsize3.value {
		return fmt.Errorf("ERROR: encode: problem w/quantization block sizes")
	}

	/* Add a Start Of Image (SOI) marker to the WSQ buffer. */
	err = e.writeShort(soi_wsq)
	if err != nil {
		return
	}
	err = e.putCNistcomWSQ(b.Dx(), b.Dy(), -1, e.quant_vals.r, o)
	if err != nil {
		return
	}
	/* Store the Wavelet filter taps to the WSQ buffer. */
	err = e.putCTransformTable(max_lofilt, max_hifilt)
	if err != nil {
		return
	}
	/* Store the quantization parameters to the WSQ buffer. */
	err = e.putCQuantizationTable()
	if err != nil {
		return
	}
	/* Store a frame header to the WSQ buffer. */
	err = e.putCFrameHeaderWSQ(b.Dx(), b.Dy(), m_shift.value, r_scale.value)
	if err != nil {
		return
	}
	/* ENCODE Block 1 */

	/* Compute Huffman table for Block 1. */

	hufftable, err = genHufftableWSQ(&e, huffbits, huffvalues, qdata, 0, []int{qsize1.value})
	if err != nil {
		return
	}

	/* Store Huffman table for Block 1 to WSQ buffer. */
	err = e.putCHuffmanTable(dht_wsq, 0, huffbits.value, huffvalues.value)
	if err != nil {
		return
	}
	/* Store Block 1's header to WSQ buffer. */
	err = e.putCBlockHeader(0)
	if err != nil {
		return
	}
	/* Compress Block 1 data. */
	err = compressBlock(&e, qdata, 0, qsize1.value, max_huffcoeff, max_huffzrun, hufftable)
	if err != nil {
		return
	}
	/* ENCODE Block 2 */

	/* Compute  Huffman table for Blocks 2 & 3. */
	hufftable, err = genHufftableWSQ(&e, huffbits, huffvalues, qdata, qsize1.value, []int{qsize2.value, qsize3.value})
	if err != nil {
		return
	}
	/* Store Huffman table for Blocks 2 & 3 to WSQ buffer. */
	err = e.putCHuffmanTable(dht_wsq, 1, huffbits.value, huffvalues.value)
	if err != nil {
		return
	}
	/* Store Block 2's header to WSQ buffer. */
	err = e.putCBlockHeader(1)
	if err != nil {
		return
	}
	/* Compress Block 2 data. */
	err = compressBlock(&e, qdata, qsize1.value, qsize2.value, max_huffcoeff, max_huffzrun, hufftable)
	if err != nil {
		return
	}
	/* ENCODE Block 3 */

	/* Store Block 3's header to WSQ buffer. */
	err = e.putCBlockHeader(1)
	if err != nil {
		return
	}
	/* Compress Block 3 data. */
	err = compressBlock(&e, qdata, qsize1.value+qsize2.value, qsize3.value, max_huffcoeff, max_huffzrun, hufftable)
	if err != nil {
		return
	}
	/* Add a End Of Image (EOI) marker to the WSQ buffer. */
	err = e.writeShort(eoi_wsq)
	if err != nil {
		return
	}

	return e.flush()
}

func wsqDecompose(encoder *encoder, fdata []float32, width, height int, hifilt []float32, hisz int, lofilt []float32, losz int) error {

	num_pix := width * height
	/* Allocate temporary floating point pixmap. */
	fdata1 := make([]float32, num_pix)

	/* Compute the Wavelet image decomposition. */
	for node := 0; node < len(encoder.wtree); node++ {
		fdataBse := (encoder.wtree[node].y * width) + encoder.wtree[node].x

		if err := getLets(fdata1, fdata, 0, fdataBse, encoder.wtree[node].leny, encoder.wtree[node].lenx,
			width, 1, hifilt, hisz, lofilt, losz, encoder.wtree[node].invrw); err != nil {
			return err
		}
		if err := getLets(fdata, fdata1, fdataBse, 0, encoder.wtree[node].lenx, encoder.wtree[node].leny,
			1, width, hifilt, hisz, lofilt, losz, encoder.wtree[node].invcl); err != nil {
			return err
		}
	}

	return nil
}

func getLets(newdata, /* image pointers for creating subband splits */
	olddata []float32,
	newIndex,
	oldIndex,
	len1, /* temporary length parameters */
	len2,
	pitch, /* pitch gives next row_col to filter */
	stride int, /*           stride gives next pixel to filter */
	hi []float32,
	hsz int, /* NEW */
	lo []float32, /* filter coefficients */
	lsz, /* NEW */
	inv int, /* spectral inversion? */
) error {
	if newdata == nil {
		return fmt.Errorf("ERROR: getLets: newdata == nil")
	}
	if olddata == nil {
		return fmt.Errorf("ERROR: getLets: olddata == nil")
	}
	if lo == nil {
		return fmt.Errorf("ERROR: getLets: lo == nil")
	}

	var lopass, hipass int /* pointers of where to put lopass
	   and hipass filter outputs */
	var p0, p1 int     /* pointers to image pixels used */
	var pix, rw_cl int /* pixel counter and row/column counter */
	var i, da_ev int   /* even or odd row/column of pixels */
	var fi_ev int
	var loc, hoc, nstr, pstr int
	var llen, hlen int
	var lpxstr, lspxstr int
	var lpx, lspx int
	var hpxstr, hspxstr int
	var hpx, hspx int
	var olle, ohle int
	var olre, ohre int
	var lle, lle2 int
	var lre, lre2 int
	var hle, hle2 int
	var hre, hre2 int

	da_ev = len2 % 2
	fi_ev = lsz % 2

	if fi_ev != 0 {
		loc = (lsz - 1) / 2
		hoc = (hsz-1)/2 - 1
		olle = 0
		ohle = 0
		olre = 0
		ohre = 0
	} else {
		loc = lsz/2 - 2
		hoc = hsz/2 - 2
		olle = 1
		ohle = 1
		olre = 1
		ohre = 1

		if loc == -1 {
			loc = 0
			olle = 0
		}
		if hoc == -1 {
			hoc = 0
			ohle = 0
		}

		for i := 0; i < hsz; i++ {
			hi[i] *= -1.0
		}
	}

	pstr = stride
	nstr = -pstr

	if da_ev != 0 {
		llen = (len2 + 1) / 2
		hlen = llen - 1
	} else {
		llen = len2 / 2
		hlen = llen
	}

	for rw_cl = 0; rw_cl < len1; rw_cl++ {
		if inv != 0 {
			hipass = newIndex + rw_cl*pitch
			lopass = hipass + hlen*stride
		} else {
			lopass = newIndex + rw_cl*pitch
			hipass = lopass + llen*stride
		}

		p0 = oldIndex + rw_cl*pitch
		p1 = p0 + (len2-1)*stride

		lspx = p0 + (loc * stride)
		lspxstr = nstr
		lle2 = olle
		lre2 = olre
		hspx = p0 + (hoc * stride)
		hspxstr = nstr
		hle2 = ohle
		hre2 = ohre
		for pix = 0; pix < hlen; pix++ {
			lpxstr = lspxstr
			lpx = lspx
			lle = lle2
			lre = lre2
			newdata[lopass] = olddata[lpx] * lo[0]
			for i = 1; i < lsz; i++ {
				if lpx == p0 {
					if lle != 0 {
						lpxstr = 0
						lle = 0
					} else {
						lpxstr = pstr
					}
				}
				if lpx == p1 {
					if lre != 0 {
						lpxstr = 0
						lre = 0
					} else {
						lpxstr = nstr
					}
				}
				lpx += lpxstr
				newdata[lopass] += olddata[lpx] * lo[i]
			}
			lopass += stride

			hpxstr = hspxstr
			hpx = hspx
			hle = hle2
			hre = hre2
			newdata[hipass] = olddata[hpx] * hi[0]
			for i = 1; i < hsz; i++ {
				if hpx == p0 {
					if hle != 0 {
						hpxstr = 0
						hle = 0
					} else {
						hpxstr = pstr
					}
				}
				if hpx == p1 {
					if hre != 0 {
						hpxstr = 0
						hre = 0
					} else {
						hpxstr = nstr
					}
				}
				hpx += hpxstr
				newdata[hipass] += olddata[hpx] * hi[i]
			}
			hipass += stride

			for i = 0; i < 2; i++ {
				if lspx == p0 {
					if lle2 != 0 {
						lspxstr = 0
						lle2 = 0
					} else {
						lspxstr = pstr
					}
				}
				lspx += lspxstr
				if hspx == p0 {
					if hle2 != 0 {
						hspxstr = 0
						hle2 = 0
					} else {
						hspxstr = pstr
					}
				}
				hspx += hspxstr
			}
		}
		if da_ev != 0 {
			lpxstr = lspxstr
			lpx = lspx
			lle = lle2
			lre = lre2
			newdata[lopass] = olddata[lpx] * lo[0]
			for i = 1; i < lsz; i++ {
				if lpx == p0 {
					if lle != 0 {
						lpxstr = 0
						lle = 0
					} else {
						lpxstr = pstr
					}
				}
				if lpx == p1 {
					if lre != 0 {
						lpxstr = 0
						lre = 0
					} else {
						lpxstr = nstr
					}
				}
				lpx += lpxstr
				newdata[lopass] += olddata[lpx] * lo[i]
			}
			lopass += stride
		}
	}
	if fi_ev == 0 {
		for i = 0; i < hsz; i++ {
			hi[i] *= -1.0
		}
	}
	return nil
}

func convertImageToFloat(m *image.Gray, m_shift, r_scale *reference[float32]) []float32 {
	bounds := m.Bounds()
	width, height := bounds.Dx(), bounds.Dy()
	var sum int
	var low, high int
	var lowDiff, highDiff float32

	fip := make([]float32, width*height)
	index := 0

	// Iterate over the image pixels
	for y := bounds.Min.Y; y < bounds.Max.Y; y++ {
		for x := bounds.Min.X; x < bounds.Max.X; x++ {
			// Get the pixel value at (x, y)
			pixelValue := int(m.GrayAt(x, y).Y)

			if pixelValue > high {
				high = pixelValue
			}
			if pixelValue < low {
				low = pixelValue
			}
			sum += pixelValue

			// Store the pixel value as a float in the fip array
			fip[index] = float32(pixelValue)
			index++
		}
	}

	m_shift.value = float32(sum) / float32(width*height)
	lowDiff = m_shift.value - float32(low)
	highDiff = float32(high) - m_shift.value

	if lowDiff >= highDiff {
		r_scale.value = lowDiff
	} else {
		r_scale.value = highDiff
	}

	r_scale.value /= 128.0

	// Scale and shift pixel values to obtain the final float array
	for i := 0; i < width*height; i++ {
		fip[i] = (fip[i] - m_shift.value) / r_scale.value
	}

	return fip
}

func variance(encoder *encoder, fip []float32, width, height int) {
	var fp int           /* temp image pointer */
	var lenx, leny int   /* dimensions of area to calculate variance */
	var skipx, skipy int /* pixels to skip to get to area for variance calculation */
	var row, col int     /* dimension counters */
	var ssq float32      /* sum of squares */
	var sum2 float32     /* variance calculation parameter */
	var sum_pix float32  /* sum of pixels */
	var vsum float32     /* variance sum for subbands 0-3 */

	vsum = 0
	for cvr := 0; cvr < 4; cvr++ {
		fp = ((encoder.qtree[cvr].y) * width) + encoder.qtree[cvr].x
		ssq = 0.0
		sum_pix = 0.0

		skipx = encoder.qtree[cvr].lenx / 8
		skipy = (9 * encoder.qtree[cvr].leny) / 32

		lenx = (3 * encoder.qtree[cvr].lenx) / 4
		leny = (7 * encoder.qtree[cvr].leny) / 16

		fp += (skipy * width) + skipx
		for row = 0; row < leny; row++ {
			for col = 0; col < lenx; col++ {
				sum_pix += fip[fp]
				ssq += fip[fp] * fip[fp]
				fp++
			}
			fp += (width - lenx)
		}
		sum2 = (sum_pix * sum_pix) / float32(lenx*leny)
		encoder.quant_vals.variance[cvr] = (ssq - sum2) / (float32(lenx*leny) - 1.0)
		vsum += encoder.quant_vals.variance[cvr]
	}

	//This part is needed to comply with WSQ 3.1
	if vsum < 20000.0 {
		for cvr := 0; cvr < num_subbands; cvr++ {
			fp = (encoder.qtree[cvr].y * width) + encoder.qtree[cvr].x
			ssq = 0
			sum_pix = 0

			lenx = encoder.qtree[cvr].lenx
			leny = encoder.qtree[cvr].leny

			for row = 0; row < leny; row++ {
				for col = 0; col < lenx; col++ {
					sum_pix += fip[fp]
					ssq += fip[fp] * fip[fp]
					fp++
				}
				fp += (width - lenx)
			}
			sum2 = (sum_pix * sum_pix) / float32(lenx*leny)
			encoder.quant_vals.variance[cvr] = (ssq - sum2) / (float32(lenx*leny) - 1.0)
		}
	} else {
		for cvr := 4; cvr < num_subbands; cvr++ {
			fp = (encoder.qtree[cvr].y * width) + encoder.qtree[cvr].x
			ssq = 0
			sum_pix = 0

			skipx = encoder.qtree[cvr].lenx / 8
			skipy = (9 * encoder.qtree[cvr].leny) / 32

			lenx = (3 * encoder.qtree[cvr].lenx) / 4
			leny = (7 * encoder.qtree[cvr].leny) / 16

			fp += (skipy * width) + skipx
			for row = 0; row < leny; row++ {
				for col = 0; col < lenx; col++ {
					sum_pix += fip[fp]
					ssq += fip[fp] * fip[fp]
					fp++
				}
				fp += (width - lenx)
			}
			sum2 = (sum_pix * sum_pix) / float32(lenx*leny)
			encoder.quant_vals.variance[cvr] = (ssq - sum2) / (float32(lenx*leny) - 1.0)
		}
	}
}

func quantize(encoder *encoder, qsize *reference[int], fip []float32, width, height int) []int {
	var row, col int                   /* temp image characteristic parameters */
	var zbin float32                   /* zero bin size */
	A := make([]float32, num_subbands) /* subband "weights" for quantization */
	m := make([]float32, num_subbands) /* subband size to image size ratios */
	/* (reciprocal of FBI spec for 'm')  */
	var m1, m2, m3 float32                 /* reciprocal constants for 'm' */
	sigma := make([]float32, num_subbands) /* square root of subband variances */
	K0 := make([]int, num_subbands)        /* initial list of subbands w/variance >= thresh */
	K1 := make([]int, num_subbands)        /* working list of subbands */
	var K, nK int                          /* pointers to sets of subbands */
	NP := make([]bool, num_subbands)       /* current subbounds with nonpositive bit rates. */
	var K0len int                          /* number of subbands in K0 */
	var Klen, nKlen int                    /* number of subbands in other subband lists */
	var NPlen int                          /* number of subbands flagged in NP */
	var S float32                          /* current frac of subbands w/positive bit rate */
	var q float32                          /* current proportionality constant */
	var P float32                          /* product of 'q/Q' ratios */

	/* Set up 'A' table. */
	for cnt := 0; cnt < strt_subband_3; cnt++ {
		A[cnt] = 1.0
	}
	A[strt_subband_3 /*52*/] = 1.32
	A[strt_subband_3+1 /*53*/] = 1.08
	A[strt_subband_3+2 /*54*/] = 1.42
	A[strt_subband_3+3 /*55*/] = 1.08
	A[strt_subband_3+4 /*56*/] = 1.32
	A[strt_subband_3+5 /*57*/] = 1.42
	A[strt_subband_3+6 /*58*/] = 1.08
	A[strt_subband_3+7 /*59*/] = 1.08

	for cnt := 0; cnt < max_subbands; cnt++ {
		encoder.quant_vals.qbss[cnt] = 0.0
		encoder.quant_vals.qzbs[cnt] = 0.0
	}

	/* Set up 'Q1' (prime) table. */
	for cnt := 0; cnt < num_subbands; cnt++ {
		if encoder.quant_vals.variance[cnt] < variance_thresh {
			encoder.quant_vals.qbss[cnt] = 0.0
		} else {
			/* NOTE: q has been taken out of the denominator in the next */
			/*       2 formulas from the original code. */
			if cnt < strt_size_region_2 /*4*/ {
				encoder.quant_vals.qbss[cnt] = 1.0
			} else {
				encoder.quant_vals.qbss[cnt] = 10.0 / (A[cnt] * float32(math.Log(float64(encoder.quant_vals.variance[cnt]))))
			}
		}
	}

	/* Set up output buffer. */
	sip := make([]int, width*height)
	var sptr int

	/* Set up 'm' table (these values are the reciprocal of 'm' in the FBI spec). */
	m1 = 1.0 / 1024.0
	m2 = 1.0 / 256.0
	m3 = 1.0 / 16.0
	for cnt := 0; cnt < strt_size_region_2; cnt++ {
		m[cnt] = m1
	}
	for cnt := strt_size_region_2; cnt < strt_size_region_3; cnt++ {
		m[cnt] = m2
	}
	for cnt := strt_size_region_3; cnt < num_subbands; cnt++ {
		m[cnt] = m3
	}

	/* Initialize 'K0' and 'K1' lists. */
	K0len = 0
	for cnt := 0; cnt < num_subbands; cnt++ {
		if encoder.quant_vals.variance[cnt] >= variance_thresh {
			K0[K0len] = cnt
			K1[K0len] = cnt
			K0len++
			/* Compute square root of subband variance. */
			sigma[cnt] = float32(math.Sqrt(float64(encoder.quant_vals.variance[cnt])))
		}
	}
	K = 0
	Klen = K0len

	for {
		/* Compute new 'S' */
		S = 0.0
		for i := 0; i < Klen; i++ {
			/* Remember 'm' is the reciprocal of spec. */
			S += m[K1[K+i]]
		}

		/* Compute product 'P' */
		P = 1.0
		for i := 0; i < Klen; i++ {
			/* Remember 'm' is the reciprocal of spec. */
			P *= float32(math.Pow(float64(sigma[K1[K+i]]/encoder.quant_vals.qbss[K1[K+i]]), float64(m[K1[K+i]])))
		}

		/* Compute new 'q' */
		q = float32((math.Pow(2, float64((encoder.quant_vals.r/S)-1.0)) / 2.5) / math.Pow(float64(P), float64(1.0/S)))

		/* Flag subbands with non-positive bitrate. */
		NP = make([]bool, num_subbands)
		NPlen = 0
		for i := 0; i < Klen; i++ {
			if (encoder.quant_vals.qbss[K1[K+i]] / q) >= (5.0 * sigma[K1[K+i]]) {
				NP[K1[K+i]] = true
				NPlen++
			}
		}

		/* If list of subbands with non-positive bitrate is empty ... */
		if NPlen == 0 {
			/* Then we are done, so break from while loop. */
			break
		}

		/* Assign new subband set to previous set K minus subbands in set NP. */
		nK = 0
		nKlen = 0
		for i := 0; i < Klen; i++ {
			if !NP[K1[K+i]] {
				K1[nK+nKlen] = K1[K+i]
				nKlen++
			}
		}

		/* Assign new set as K. */
		K = nK
		Klen = nKlen
	}

	/* Flag subbands that are in set 'K0' (the very first set). */
	nK = 0

	for i := 0; i < K0len; i++ {
		K1[nK+K0[i]] = 1 /* MO: was = TRUE */
	}
	/* Set 'Q' values. */

	for cnt := 0; cnt < num_subbands; cnt++ {
		if K1[nK+cnt] != 0 {
			encoder.quant_vals.qbss[cnt] /= q
		} else {
			encoder.quant_vals.qbss[cnt] = 0.0
		}
		encoder.quant_vals.qzbs[cnt] = 1.2 * encoder.quant_vals.qbss[cnt]
	}

	/* Now ready to compute and store bin widths for subbands. */
	for cnt := 0; cnt < num_subbands; cnt++ {
		fptr := (encoder.qtree[cnt].y * width) + encoder.qtree[cnt].x

		if encoder.quant_vals.qbss[cnt] != 0.0 {

			zbin = encoder.quant_vals.qzbs[cnt] / 2.0

			for row = 0; row < encoder.qtree[cnt].leny; row++ {
				for col = 0; col < encoder.qtree[cnt].lenx; col++ {
					if -zbin <= fip[fptr] && fip[fptr] <= zbin {
						sip[sptr] = 0
					} else if fip[fptr] > 0.0 {
						sip[sptr] = int(((fip[fptr] - zbin) / encoder.quant_vals.qbss[cnt]) + 1.0)
					} else {
						sip[sptr] = int(((fip[fptr] + zbin) / encoder.quant_vals.qbss[cnt]) - 1.0)
					}
					sptr++
					fptr++
				}
				fptr += width - encoder.qtree[cnt].lenx
			}
		}
	}

	qsize.value = sptr

	return sip
}

func quantBlockSizes(encoder *encoder, oqsize1, oqsize2, oqsize3 *reference[int]) {
	var qsize1, qsize2, qsize3 int
	var node int

	/* Compute temporary sizes of 3 WSQ subband blocks. */
	qsize1 = encoder.wtree[14].lenx * encoder.wtree[14].leny
	qsize2 = (encoder.wtree[5].leny * encoder.wtree[1].lenx) +
		(encoder.wtree[4].lenx * encoder.wtree[4].leny)
	qsize3 = (encoder.wtree[2].lenx * encoder.wtree[2].leny) +
		(encoder.wtree[3].lenx * encoder.wtree[3].leny)

	/* Adjust size of quantized WSQ subband blocks. */
	for node = 0; node < strt_size_region_2; node++ {
		if encoder.quant_vals.qbss[node] == 0.0 {
			qsize1 -= (encoder.qtree[node].lenx * encoder.qtree[node].leny)
		}
	}

	for node = strt_subband_2; node < strt_subband_3; node++ {
		if encoder.quant_vals.qbss[node] == 0.0 {
			qsize2 -= (encoder.qtree[node].lenx * encoder.qtree[node].leny)
		}
	}

	for node = strt_subband_2; node < strt_subband_del; node++ {
		if encoder.quant_vals.qbss[node] == 0.0 {
			qsize3 -= (encoder.qtree[node].lenx * encoder.qtree[node].leny)
		}
	}

	oqsize1.value = qsize1
	oqsize2.value = qsize2
	oqsize3.value = qsize3
}

func genHufftableWSQ(encoder *encoder, ohuffbits, ohuffvalues *reference[[]int], sip []int, offset int, block_sizes []int) (hufftable2 []huffCode, err error) {
	var codesize []int        /* code sizes to use */
	var tempSize int          /* last huffvalue */
	var huffbits []int        /* huffbits values */
	var huffvalues []int      /* huffvalues */
	var huffcounts []int      /* counts for each huffman category */
	var huffcounts2 []int     /* counts for each huffman category */
	var hufftable1 []huffCode /* hufftables */

	huffcounts, err = countBlock(max_huffcounts_wsq, sip, offset, block_sizes[0], max_huffcoeff, max_huffzrun)
	if err != nil {
		return
	}

	for i := 1; i < len(block_sizes); i++ {
		huffcounts2, err = countBlock(max_huffcounts_wsq, sip, offset+block_sizes[i-1], block_sizes[i], max_huffcoeff, max_huffzrun)
		if err != nil {
			return
		}

		for j := 0; j < max_huffcounts_wsq; j++ {
			huffcounts[j] += huffcounts2[j]
		}
	}

	codesize = findHuffSizes(huffcounts, max_huffcounts_wsq)

	/* tells if codesize is greater than MAX_HUFFBITS */
	adjust := &reference[bool]{false}

	huffbits = findNumHuffSizes(adjust, codesize, max_huffcounts_wsq)

	if adjust.value {
		if err := sortHuffbits(huffbits); err != nil {
			return nil, err
		}
	}

	huffvalues = sortCodeSizes(codesize, max_huffcounts_wsq)

	hufftable1, tempSize = buildHuffsizes(huffbits, max_huffcounts_wsq)
	buildHuffcodes(hufftable1)
	checkHuffcodesWSQ(hufftable1, tempSize)

	hufftable2 = buildHuffcodeTable(hufftable1, tempSize, huffvalues, max_huffcounts_wsq)

	ohuffbits.value = huffbits
	ohuffvalues.value = huffvalues
	return hufftable2, nil
}

func countBlock(
	// int **ocounts,     /* output count for each huffman catetory */
	max_huffcounts int, /* maximum number of counts */
	sip []int, /* quantized data */
	sip_offset int, /* offset into sip */
	sip_siz int, /* size of block being compressed */
	MaxCoeff int, /* maximum values for coefficients */
	MaxZRun int, /* maximum zero runs */
) ([]int, error) {
	var counts []int    /* count for each huffman category */
	var LoMaxCoeff int  /* lower (negative) MaxCoeff limit */
	var pix int         /* temp pixel pointer */
	var rcnt, state int /* zero run count and if current pixel
	is in a zero run or just a coefficient */
	var cnt int /* pixel counter */

	if MaxCoeff < 0 || MaxCoeff > 0xffff {
		return nil, fmt.Errorf("ERROR : countBlock : MaxCoeff out of range.")
	}
	if MaxZRun < 0 || MaxZRun > 0xffff {
		return nil, fmt.Errorf("ERROR : countBlock : MaxZRun out of range.")
	}
	/* Ininitalize vector of counts to 0. */
	counts = make([]int, max_huffcounts+1)
	/* Set last count to 1. */
	counts[max_huffcounts] = 1

	LoMaxCoeff = 1 - MaxCoeff
	state = coeff_code
	for cnt = sip_offset; cnt < sip_siz; cnt++ {
		pix = sip[cnt]
		switch state {

		case coeff_code: /* for runs of zeros */
			if pix == 0 {
				state = run_code
				rcnt = 1
				break
			}
			if pix > MaxCoeff {
				if pix > 255 {
					counts[103]++ /* 16bit pos esc */
				} else {
					counts[101]++ /* 8bit pos esc */
				}
			} else if pix < LoMaxCoeff {
				if pix < -255 {
					counts[104]++ /* 16bit neg esc */
				} else {
					counts[102]++ /* 8bit neg esc */
				}
			} else {
				counts[pix+180]++ /* within table */
			}
			break

		case run_code: /* get length of zero run */
			if pix == 0 && rcnt < 0xFFFF {
				rcnt++
				break
			}
			/* limit rcnt to avoid EOF problem in bitio.c */
			if rcnt <= MaxZRun {
				counts[rcnt]++ /** log zero run length **/
			} else if rcnt <= 0xFF {
				counts[105]++
			} else if rcnt <= 0xFFFF {
				counts[106]++ /* 16bit zrun esc */
			} else {
				return nil, fmt.Errorf("ERROR: countBlock : Zrun to long in count block.")
			}

			if pix != 0 {
				if pix > MaxCoeff { /** log current pix **/
					if pix > 255 {
						counts[103]++ /* 16bit pos esc */
					} else {
						counts[101]++ /* 8bit pos esc */
					}
				} else if pix < LoMaxCoeff {
					if pix < -255 {
						counts[104]++ /* 16bit neg esc */
					} else {
						counts[102]++ /* 8bit neg esc */
					}
				} else {
					counts[pix+180]++ /* within table */
				}
				state = coeff_code
			} else {
				rcnt = 1
				state = run_code
			}
			break
		}
	}
	if state == run_code { /** log zero run length **/
		if rcnt <= MaxZRun {
			counts[rcnt]++
		} else if rcnt <= 0xFF {
			counts[105]++
		} else if rcnt <= 0xFFFF {
			counts[106]++ /* 16bit zrun esc */
		} else {
			return nil, fmt.Errorf("ERROR: countBlock : Zrun to long in count block.")
		}
	}

	return counts, nil
}

func findNumHuffSizes(adjust *reference[bool], codesize []int, max_huffcounts int) []int {
	adjust.value = false

	/* Allocate 2X desired number of bits due to possible codesize. */
	bits := make([]int, 2*max_huffbits)

	for i := 0; i < max_huffcounts; i++ {
		if codesize[i] != 0 {
			bits[codesize[i]-1]++
		}
		if codesize[i] > max_huffbits {
			adjust.value = true
		}
	}
	return bits
}

func sortCodeSizes(codesize []int, max_huffcounts int) []int {
	/* defines order of huffman codelengths in relation to the code sizes */
	values := make([]int, max_huffcounts+1)
	var i2 int
	for i := 1; i <= (max_huffbits << 1); i++ {
		for i3 := 0; i3 < max_huffcounts; i3++ {
			if codesize[i3] == i {
				values[i2] = i3
				i2++
			}
		}
	}
	return values
}

func checkHuffcodesWSQ(hufftable []huffCode, last_size int) {
	var all_ones bool

	for i := 0; i < last_size; i++ {
		all_ones = true
		for k := 0; (k < hufftable[i].size) && all_ones; k++ {
			all_ones = (all_ones && (((hufftable[i].code >> k) & 0x0001) != 0))
		}
		if all_ones {
			log.Println("WARNING: checkHuffcodesWSQ: A code in the hufftable contains an all 1's code. This image may still be decodable. It is not compliant with the WSQ specification.")
		}
	}

}

func buildHuffcodeTable(in_huffcode_table []huffCode, last_size int, values []int, max_huffcounts int) []huffCode {
	new_huffcode_table := make([]huffCode, max_huffcounts+1)
	for i := 0; i < len(new_huffcode_table); i++ {
		new_huffcode_table[i] = huffCode{}
	}

	for size := 0; size < last_size; size++ {
		new_huffcode_table[values[size]].code = in_huffcode_table[size].code
		new_huffcode_table[values[size]].size = in_huffcode_table[size].size
	}

	return new_huffcode_table
}

func sortHuffbits(bits []int) error {
	var i, j int
	var l1, l2, l3 int

	l3 = max_huffbits << 1 /* 32 */
	l1 = l3 - 1            /* 31 */
	l2 = max_huffbits - 1  /* 15 */

	tbits := make([]int, l3)

	for i = 0; i < max_huffbits<<1; i++ {
		tbits[i] = bits[i]
	}

	for i = l1; i > l2; i-- {
		for tbits[i] > 0 {
			j = i - 2
			for tbits[j] == 0 {
				j--
			}
			tbits[i] -= 2
			tbits[i-1] += 1
			tbits[j+1] += 2
			tbits[j] -= 1
		}
		tbits[i] = 0
	}

	for tbits[i] == 0 {
		i--
	}

	tbits[i] -= 1

	for i = 0; i < max_huffbits<<1; i++ {
		bits[i] = tbits[i] // int(tbits[i] & 0xff)
	}

	for i = max_huffbits; i < l3; i++ {
		if bits[i] > 0 {
			return fmt.Errorf("ERROR: sortHuffbits: Code length is greater than 16.")
		}
	}

	return nil
}

func findHuffSizes(freq []int, max_huffcounts int) []int {
	var value1 int
	/* smallest and next smallest frequency */
	var value2 int /* of difference occurrence in the largest difference category */

	/* codesizes for each category */
	codesize := make([]int, max_huffcounts+1)

	/* pointer used to generate codesizes */
	others := make([]int, max_huffcounts+1)

	for i := 0; i <= max_huffcounts; i++ {
		others[i] = -1
	}

	for {

		values := findLeastFreq(freq, max_huffcounts)
		value1 = values[0]
		value2 = values[1]

		if value2 == -1 {
			break
		}

		freq[value1] += freq[value2]
		freq[value2] = 0

		codesize[value1]++
		for others[value1] != -1 {
			value1 = others[value1]
			codesize[value1]++
		}
		others[value1] = value2
		codesize[value2]++

		for others[value2] != -1 {
			value2 = others[value2]
			codesize[value2]++
		}
	}

	return codesize
}

func findLeastFreq(freq []int, max_huffcounts int) []int {
	var code_temp int           /*store code*/
	var value_temp int          /*store size*/
	code2 := int(^uint(0) >> 1) /*next smallest frequency in largest diff category*/
	code1 := int(^uint(0) >> 1) /*smallest frequency in largest difference category*/
	set := 1                    /*flag first two non-zero frequency values*/

	value1 := -1
	value2 := -1

	for i := 0; i <= max_huffcounts; i++ {
		if freq[i] == 0 {
			continue
		}
		if set == 1 {
			code1 = freq[i]
			value1 = i
			set++
			continue
		}
		if set == 2 {
			code2 = freq[i]
			value2 = i
			set++
		}
		code_temp = freq[i]
		value_temp = i
		if code1 < code_temp && code2 < code_temp {
			continue
		}
		if (code_temp < code1) || (code_temp == code1 && value_temp > value1) {
			code2 = code1
			value2 = value1
			code1 = code_temp
			value1 = value_temp
			continue
		}
		if (code_temp < code2) || (code_temp == code2 && value_temp > value2) {
			code2 = code_temp
			value2 = value_temp
		}
	}
	return []int{value1, value2}
}

func compressBlock(encoder *encoder,
	sip []int, /* quantized image */
	offset,
	length,
	MaxCoeff, /* Maximum values for coefficients  */
	MaxZRun int, /* Maximum zero runs */
	codes []huffCode, /* huffman code table  */
) error {
	var LoMaxCoeff int  /* lower (negative) MaxCoeff limit */
	var pix int         /* temp pixel pointer */
	var rcnt, state int /* zero run count and if current pixel
	is in a zero run or just a coefficient */
	var cnt int /* pixel counter */

	if MaxCoeff < 0 || MaxCoeff > 0xffff {
		return fmt.Errorf("ERROR: compressBlock: MaxCoeff out of range.")
	}
	if MaxZRun < 0 || MaxZRun > 0xffff {
		return fmt.Errorf("ERROR: compressBlock: MaxZRun out of range.")
	}
	LoMaxCoeff = 1 - MaxCoeff

	outbit := &reference[int]{7}
	bytes := &reference[int]{0}
	bits := &reference[int]{0}

	state = coeff_code
	for cnt = offset; cnt < length; cnt++ {
		pix = sip[cnt]

		switch state {

		case coeff_code:
			if pix == 0 {
				state = run_code
				rcnt = 1
				break
			}
			if pix > MaxCoeff {
				if pix > 255 {
					/* 16bit pos esc */
					if err := writeBits(encoder, codes[103].size, codes[103].code, outbit, bits, bytes); err != nil {
						return err
					}
					if err := writeBits(encoder, 16, pix, outbit, bits, bytes); err != nil {
						return err
					}
				} else {
					/* 8bit pos esc */
					if err := writeBits(encoder, codes[101].size, codes[101].code, outbit, bits, bytes); err != nil {
						return err
					}
					if err := writeBits(encoder, 8, pix, outbit, bits, bytes); err != nil {
						return err
					}
				}
			} else if pix < LoMaxCoeff {
				if pix < -255 {
					/* 16bit neg esc */
					if err := writeBits(encoder, codes[104].size, codes[104].code, outbit, bits, bytes); err != nil {
						return err
					}
					if err := writeBits(encoder, 16, -(pix), outbit, bits, bytes); err != nil {
						return err
					}
				} else {
					/* 8bit neg esc */
					if err := writeBits(encoder, codes[102].size, codes[102].code, outbit, bits, bytes); err != nil {
						return err
					}
					if err := writeBits(encoder, 8, -(pix), outbit, bits, bytes); err != nil {
						return err
					}
				}
			} else {
				/* within table */
				if err := writeBits(encoder, codes[pix+180].size, codes[pix+180].code, outbit, bits, bytes); err != nil {
					return err
				}
			}
			break

		case run_code:
			if pix == 0 && rcnt < 0xFFFF {
				rcnt++
				break
			}
			if rcnt <= MaxZRun {
				/* log zero run length */
				if err := writeBits(encoder, codes[rcnt].size, codes[rcnt].code, outbit, bits, bytes); err != nil {
					return err
				}
			} else if rcnt <= 0xFF {
				/* 8bit zrun esc */
				if err := writeBits(encoder, codes[105].size, codes[105].code, outbit, bits, bytes); err != nil {
					return err
				}
				if err := writeBits(encoder, 8, rcnt, outbit, bits, bytes); err != nil {
					return err
				}
			} else if rcnt <= 0xFFFF {
				/* 16bit zrun esc */
				if err := writeBits(encoder, codes[106].size, codes[106].code, outbit, bits, bytes); err != nil {
					return err
				}
				if err := writeBits(encoder, 16, rcnt, outbit, bits, bytes); err != nil {
					return err
				}
			} else {
				return fmt.Errorf("ERROR : compress_block : zrun too large.")
			}

			if pix != 0 {
				if pix > MaxCoeff {
					/** log current pix **/
					if pix > 255 {
						/* 16bit pos esc */
						if err := writeBits(encoder, codes[103].size, codes[103].code, outbit, bits, bytes); err != nil {
							return err
						}
						if err := writeBits(encoder, 16, pix, outbit, bits, bytes); err != nil {
							return err
						}
					} else {
						/* 8bit pos esc */
						if err := writeBits(encoder, codes[101].size, codes[101].code, outbit, bits, bytes); err != nil {
							return err
						}
						if err := writeBits(encoder, 8, pix, outbit, bits, bytes); err != nil {
							return err
						}
					}
				} else if pix < LoMaxCoeff {
					if pix < -255 {
						/* 16bit neg esc */
						if err := writeBits(encoder, codes[104].size, codes[104].code, outbit, bits, bytes); err != nil {
							return err
						}
						if err := writeBits(encoder, 16, -pix, outbit, bits, bytes); err != nil {
							return err
						}
					} else {
						/* 8bit neg esc */
						if err := writeBits(encoder, codes[102].size, codes[102].code, outbit, bits, bytes); err != nil {
							return err
						}
						if err := writeBits(encoder, 8, -pix, outbit, bits, bytes); err != nil {
							return err
						}
					}
				} else {
					/* within table */
					if err := writeBits(encoder, codes[pix+180].size, codes[pix+180].code, outbit, bits, bytes); err != nil {
						return err
					}
				}
				state = coeff_code
			} else {
				rcnt = 1
				state = run_code
			}
			break
		}
	}
	if state == run_code {
		if rcnt <= MaxZRun {
			if err := writeBits(encoder, codes[rcnt].size, codes[rcnt].code, outbit, bits, bytes); err != nil {
				return err
			}
		} else if rcnt <= 0xFF {
			if err := writeBits(encoder, codes[105].size, codes[105].code, outbit, bits, bytes); err != nil {
				return err
			}
			if err := writeBits(encoder, 8, rcnt, outbit, bits, bytes); err != nil {
				return err
			}
		} else if rcnt <= 0xFFFF {
			if err := writeBits(encoder, codes[106].size, codes[106].code, outbit, bits, bytes); err != nil {
				return err
			}
			if err := writeBits(encoder, 16, rcnt, outbit, bits, bytes); err != nil {
				return err
			}
		} else {
			return fmt.Errorf("ERROR: compressBlock: zrun2 too large.")
		}
	}

	return flushBits(encoder, outbit, bits, bytes)
}

func writeBits(
	encoder *encoder,
	size, /* numbers bits of code to write into buffer   */
	code int, /* info to write into buffer                   */
	outbit, /* current bit location in out buffer byte     */
	bits, /* byte to write to output buffer              */
	bytes *reference[int], /* count of number bytes written to the buffer */
) error {
	num := size
	for num--; num >= 0; num-- {
		bits.value <<= 1
		bits.value |= (code >> num) & 0x0001

		if outbit.value-1 < 0 {
			if err := encoder.write([]byte{byte(bits.value)}); err != nil {
				return err
			}
			if (bits.value & 0xFF) == 0xFF {
				if err := encoder.write([]byte{0}); err != nil {
					return err
				}
				bytes.value++
			}
			bytes.value++
			outbit.value = 7
			bits.value = 0
		} else {
			outbit.value--
		}
	}
	return nil
}

func flushBits(
	encoder *encoder, /* output data buffer */
	outbit, /* current bit location in out buffer byte */
	bits, /* byte to write to output buffer */
	bytes *reference[int], /* count of number bytes written to the buffer */
) error {
	var cnt int /* temp counter */

	if outbit.value != 7 {
		for cnt >= 0 {
			bits.value <<= 1
			bits.value |= 0x01
			cnt--
		}

		if err := encoder.write([]byte{byte(bits.value)}); err != nil {
			return err
		}
		if bits.value == 0xFF {
			bits.value = 0
			encoder.write([]byte{0})
			bytes.value++
		}
		bytes.value++
		outbit.value = 7
		bits.value = 0
	}
	return nil
}

func fetToString(fet map[string]string) string {
	var result strings.Builder

	for key, value := range fet {
		if key == "" || value == "" {
			continue
		}

		encodedKey := url.QueryEscape(key)
		encodedValue := url.QueryEscape(value)

		result.WriteString(encodedKey)
		result.WriteString(" ")
		result.WriteString(encodedValue)
		result.WriteString("\n")
	}

	return result.String()
}
