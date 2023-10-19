package wsq

import (
	"fmt"
	"io"
	"math"
)

const (
	NCM_EXT           = "ncm"
	NCM_HEADER        = "NIST_COM"       /* mandatory */
	NCM_PIX_WIDTH     = "PIX_WIDTH"      /* mandatory */
	NCM_PIX_HEIGHT    = "PIX_HEIGHT"     /* mandatory */
	NCM_PIX_DEPTH     = "PIX_DEPTH"      /* 1,8,24 (mandatory)*/
	NCM_PPI           = "PPI"            /* -1 if unknown (mandatory)*/
	NCM_COLORSPACE    = "COLORSPACE"     /* RGB,YCbCr,GRAY */
	NCM_N_CMPNTS      = "NUM_COMPONENTS" /* [1..4] (mandatory w/hv_factors)*/
	NCM_HV_FCTRS      = "HV_FACTORS"     /* H0,V0:H1,V1:...*/
	NCM_INTRLV        = "INTERLEAVE"     /* 0,1 (mandatory w/depth=24) */
	NCM_COMPRESSION   = "COMPRESSION"    /* NONE,JPEGB,JPEGL,WSQ */
	NCM_JPEGB_QUAL    = "JPEGB_QUALITY"  /* [20..95] */
	NCM_JPEGL_PREDICT = "JPEGL_PREDICT"  /* [1..7] */
	NCM_WSQ_RATE      = "WSQ_BITRATE"    /* ex. .75,2.25 (-1.0 if unknown)*/
	NCM_LOSSY         = "LOSSY"          /* 0,1 */

	NCM_HISTORY    = "HISTORY"    /* ex. SD historical data */
	NCM_FING_CLASS = "FING_CLASS" /* ex. A,L,R,S,T,W */
	NCM_SEX        = "SEX"        /* m,f */
	NCM_SCAN_TYPE  = "SCAN_TYPE"  /* l,i */
	NCM_FACE_POS   = "FACE_POS"   /* f,p */
	NCM_AGE        = "AGE"
	NCM_SD_ID      = "SD_ID" /* 4,9,10,14,18 */
)

type writer interface {
	Flush() error
	io.Writer
	io.ByteWriter
}

type encoder struct {
	token
	w          writer
	quant_vals quantization
}

func (e *encoder) init() {
	e.quant_vals = quantization{}
	e.tableDHT = make([]*tableDHT, max_dht_tables)
	for i := 0; i < max_dht_tables; i++ {
		e.tableDHT[i] = new(tableDHT)
		e.tableDHT[i].tabdef = 0
		e.tableDHT[i].huffbits = make([]int, max_huffbits)
		e.tableDHT[i].huffvalues = make([]int, max_huffcounts_wsq+1)
	}
	e.tableDTT.lofilt = lo_filt_not_even_8x8_1
	e.tableDTT.hifilt = hi_filt_not_even_8x8_1
}

func (e *encoder) putCNistcomWSQ(width, height, ppi int, bitrate float32, o *Options) error {
	if o == nil {
		return nil
	}
	nistcom := make(map[string]string)
	//These attributes will be filled later
	nistcom[NCM_HEADER] = "---"
	nistcom[NCM_PIX_WIDTH] = "---"
	nistcom[NCM_PIX_HEIGHT] = "---"
	nistcom[NCM_PIX_DEPTH] = "---"
	nistcom[NCM_PPI] = "---"
	nistcom[NCM_LOSSY] = "---"
	nistcom[NCM_COLORSPACE] = "---"
	nistcom[NCM_COMPRESSION] = "---"
	nistcom[NCM_WSQ_RATE] = "---"

	if o.Metadata != nil {
		for k, v := range o.Metadata {
			nistcom[k] = v
		}
	}

	nistcom[NCM_HEADER] = fmt.Sprintf("%d", len(nistcom))
	nistcom[NCM_PIX_WIDTH] = fmt.Sprintf("%d", width)
	nistcom[NCM_PIX_HEIGHT] = fmt.Sprintf("%d", height)
	nistcom[NCM_PPI] = fmt.Sprintf("%d", ppi)
	nistcom[NCM_PIX_DEPTH] = "8" //WSQ has always 8 bpp
	nistcom[NCM_LOSSY] = "1"     //WSQ is always lossy
	nistcom[NCM_COLORSPACE] = "GRAY"
	nistcom[NCM_COMPRESSION] = "WSQ"
	nistcom[NCM_WSQ_RATE] = fmt.Sprintf("%f", bitrate)

	if err := e.putCComment(com_wsq, fetToString(nistcom)); err != nil {
		return err
	}
	if o.Comments != nil {
		for _, s := range o.Comments {
			if s != "" {
				if err := e.putCComment(com_wsq, s); err != nil {
					return err
				}
			}
		}
	}

	return nil
}

func (e *encoder) putCComment(marker int, comment string) error {
	if err := e.writeShort(marker); err != nil {
		return err
	}

	/* comment size */
	hdr_size := 2 + len(comment)
	if err := e.writeShort(hdr_size); err != nil {
		return err
	}

	return e.write([]byte(comment)) // FIXME: should be UTF-8. Check FBI spec.
}

func (e *encoder) putCTransformTable(losz, hisz int) error {
	var coef int64         /* filter coefficient indicator */
	var int_dat int64      /* temp variable */
	var dbl_tmp float32    /* temp variable */
	var scale_ex, sign int /* exponent scaling and sign parameters */

	/* FIXME how big can hisz,losz get */
	if losz < 0 || losz > int(^uint(0)>>1)/2 {
		return fmt.Errorf("ERROR: putCTransformTable: Writing transform table: losz out of range")
	}
	if hisz < 0 || hisz > int(^uint(0)>>1)/2 {
		return fmt.Errorf("ERROR: putCTransformTable: Writing transform table: hisz out of range")
	}
	if err := e.writeShort(dtt_wsq); err != nil {
		return err
	}

	/* table size */
	if err := e.writeShort(58); err != nil {
		return err
	}

	/* number analysis lowpass coefficients */
	if err := e.writeByte(losz); err != nil {
		return err
	}

	/* number analysis highpass coefficients */
	if err := e.writeByte(hisz); err != nil {
		return err
	}

	for coef = int64(losz >> 1); float64(coef) < float64(losz); coef++ {
		dbl_tmp = e.tableDTT.lofilt[coef]
		if dbl_tmp >= 0.0 {
			sign = 0
		} else {
			sign = 1
			dbl_tmp *= -1.0
		}
		scale_ex = 0
		if dbl_tmp == 0.0 {
			int_dat = 0
		} else if dbl_tmp < 4294967295.0 {
			for dbl_tmp < 4294967295.0 {
				scale_ex += 1
				dbl_tmp *= 10.0
			}
			scale_ex -= 1
			int_dat = int64(math.Round(float64(dbl_tmp / 10.0)))
		} else {
			dbl_tmp = e.tableDTT.lofilt[coef]
			return fmt.Errorf("ERROR: putCTransformTable: lofilt[%d] to high at %f", coef, dbl_tmp)
		}

		if err := e.writeByte(sign); err != nil {
			return err
		}
		if err := e.writeByte(scale_ex); err != nil {
			return err
		}
		if err := e.writeLong(int_dat); err != nil {
			return err
		}
	}

	for coef = int64(hisz >> 1); int64(coef&0xFFFFFFFF) < int64(hisz); coef++ {
		dbl_tmp = e.tableDTT.hifilt[coef]
		if dbl_tmp >= 0.0 {
			sign = 0
		} else {
			sign = 1
			dbl_tmp *= -1.0
		}
		scale_ex = 0
		if dbl_tmp == 0.0 {
			int_dat = 0
		} else if dbl_tmp < 4294967295.0 {
			for dbl_tmp < 4294967295.0 {
				scale_ex += 1
				dbl_tmp *= 10.0
			}
			scale_ex -= 1
			int_dat = int64(math.Round(float64(dbl_tmp / 10.0)))
		} else {
			dbl_tmp = e.tableDTT.hifilt[coef]
			return fmt.Errorf("ERROR: putCTransformTable: hifilt[%d] to high at %f", coef, dbl_tmp)
		}
		if err := e.writeByte(sign); err != nil {
			return err
		}
		if err := e.writeByte(scale_ex); err != nil {
			return err
		}
		if err := e.writeLong(int_dat); err != nil {
			return err
		}
	}
	return nil
}
func (e *encoder) putCQuantizationTable() error {
	var scale_ex, scale_ex2 int /* exponent scaling parameters */
	var shrt_dat, shrt_dat2 int /* temp variables */
	var flt_tmp float32         /* temp variable */

	if err := e.writeShort(dqt_wsq); err != nil {
		return err
	}

	/* table size */
	if err := e.writeShort(389); err != nil {
		return err
	}

	/* exponent scaling value */
	if err := e.writeByte(2); err != nil {
		return err
	}

	/* quantizer bin center parameter */
	if err := e.writeShort(44); err != nil {
		return err
	}

	for sub := 0; sub < 64; sub++ {
		if sub >= 0 && sub < 60 {
			if e.quant_vals.qbss[sub] != 0.0 {
				flt_tmp = e.quant_vals.qbss[sub]
				scale_ex = 0
				if flt_tmp < 65535 {
					for flt_tmp < 65535 {
						scale_ex += 1
						flt_tmp *= 10
					}
					scale_ex -= 1
					shrt_dat = int(math.Round(float64(flt_tmp) / 10.0))
				} else {
					flt_tmp = e.quant_vals.qbss[sub]
					return fmt.Errorf("ERROR: putc_quantization_table: Q[%d] to high at %f", sub, flt_tmp)
				}

				flt_tmp = e.quant_vals.qzbs[sub]
				scale_ex2 = 0
				if flt_tmp < 65535 {
					for flt_tmp < 65535 {
						scale_ex2 += 1
						flt_tmp *= 10
					}
					scale_ex2 -= 1
					shrt_dat2 = int(math.Round(float64(flt_tmp) / 10.0))
				} else {
					flt_tmp = e.quant_vals.qzbs[sub]

					return fmt.Errorf("ERROR: putc_quantization_table: Z[%d] to high at %f", sub, flt_tmp)
				}
			} else {
				scale_ex = 0
				scale_ex2 = 0
				shrt_dat = 0
				shrt_dat2 = 0
			}
		} else {
			scale_ex = 0
			scale_ex2 = 0
			shrt_dat = 0
			shrt_dat2 = 0
		}

		if err := e.writeByte(scale_ex); err != nil {
			return err
		}
		if err := e.writeShort(shrt_dat); err != nil {
			return err
		}
		if err := e.writeByte(scale_ex2); err != nil {
			return err
		}
		if err := e.writeShort(shrt_dat2); err != nil {
			return err
		}
	}

	return nil
}

func (e *encoder) putCFrameHeaderWSQ(width, height int, m_shift, r_scale float32) error {
	var flt_tmp float32 /* temp variable */
	var scale_ex int    /* exponent scaling parameter */

	var shrt_dat int /* temp variable */
	/* +2 = 2 */
	if err := e.writeShort(sof_wsq); err != nil {
		return err
	}

	/* size of frame header */
	if err := e.writeShort(17); err != nil {
		return err
	}

	/* black pixel */
	/* +1 = 3 */
	if err := e.writeByte(0); err != nil {
		return err
	}

	/* white pixel */
	/* +1 = 4 */
	if err := e.writeByte(255); err != nil {
		return err
	}
	/* +2 = 5 */
	if err := e.writeShort(height); err != nil {
		return err
	}
	/* +2 = 7 */
	if err := e.writeShort(width); err != nil {
		return err
	}

	flt_tmp = m_shift
	scale_ex = 0
	if flt_tmp != 0.0 {
		for flt_tmp < 65535 {
			scale_ex += 1
			flt_tmp *= 10
		}
		scale_ex -= 1
		shrt_dat = int(math.Round(float64(flt_tmp) / 10.0))
	} else {
		shrt_dat = 0
	}
	/* +1 = 9 */
	if err := e.writeByte(scale_ex); err != nil {
		return err
	}
	/* +2 = 11 */
	if err := e.writeShort(shrt_dat); err != nil {
		return err
	}

	flt_tmp = r_scale
	scale_ex = 0
	if flt_tmp != 0.0 {
		for flt_tmp < 65535 {
			scale_ex += 1
			flt_tmp *= 10
		}
		scale_ex -= 1
		shrt_dat = int(math.Round(float64(flt_tmp) / 10.0))
	} else {
		shrt_dat = 0
	}
	/* +1 = 12 */
	if err := e.writeByte(scale_ex); err != nil {
		return err
	}
	/* +2 = 13 */
	if err := e.writeShort(shrt_dat); err != nil {
		return err
	}
	/* +1 = 15 */
	if err := e.writeByte(0); err != nil {
		return err
	}
	/* +2 = 17 */
	if err := e.writeShort(0); err != nil {
		return err
	}
	return nil
}

func (e *encoder) putCHuffmanTable(marker, tableId int, huffbits, huffvalues []int) error {
	if err := e.writeShort(marker); err != nil {
		return err
	}

	/* "value(2) + table id(1) + bits(16)" */
	table_len := 3 + max_huffbits
	values_offset := table_len
	for i := 0; i < max_huffbits; i++ {
		table_len += huffbits[i] /* values size */
	}

	/* Table Len */
	if err := e.writeShort(table_len); err != nil {
		return err
	}

	/* Table ID */
	if err := e.writeByte(tableId); err != nil {
		return err
	}

	/* Huffbits (MAX_HUFFBITS) */
	for i := 0; i < max_huffbits; i++ {
		if err := e.writeByte(huffbits[i]); err != nil {
			return err
		}
	}

	/* Huffvalues (MAX_HUFFCOUNTS) */
	for i := 0; i < table_len-values_offset; i++ {
		if err := e.writeByte(huffvalues[i]); err != nil {
			return err
		}
	}

	return nil
}
func (e *encoder) putCBlockHeader(table int) error {
	if err := e.writeShort(sob_wsq); err != nil {
		return err
	}

	/* block header size */
	if err := e.writeShort(3); err != nil {
		return err
	}

	return e.writeByte(table)
}

func (e *encoder) write(p []byte) error {
	_, err := e.w.Write(p)
	if err != nil {
		return err
	}
	return nil
}

func (e *encoder) writeByte(value int) error {
	return e.w.WriteByte(byte(value & 0xFF))
}

func (e *encoder) writeShort(value int) error {
	b1 := byte((value >> 8) & 0xFF)
	b2 := byte(value & 0xFF)
	return e.write([]byte{b1, b2})
}

func (e *encoder) writeLong(value int64) error {
	b1 := byte((value >> 24) & 0xFF)
	b2 := byte((value >> 16) & 0xFF)
	b3 := byte((value >> 8) & 0xFF)
	b4 := byte(value & 0xFF)

	err := e.write([]byte{b1, b2, b3, b4})
	if err != nil {
		return err
	}

	return nil
}

func (e *encoder) flush() error {
	return e.w.Flush()
}
