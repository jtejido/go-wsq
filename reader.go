package wsq

import (
	"fmt"
	"image"
	"image/color"
	"io"
)

func Decode(r io.Reader) (image.Image, error) {
	return decode(r)
}

func decode(r io.Reader) (image.Image, error) {
	var d decoder
	d.r = r
	d.init()
	/* Read the SOI marker. */
	_, err := d.getCMarkerWSQ(soi_wsq)
	if err != nil {
		return nil, err
	}

	/* Read in supporting tables up to the SOF marker. */
	marker, err := d.getCMarkerWSQ(tbls_n_sof)
	if err != nil {
		return nil, err
	}
	for marker != sof_wsq {
		if err := d.getCTableWSQ(marker); err != nil {
			return nil, err
		}
		marker, err = d.getCMarkerWSQ(tbls_n_sof)
		if err != nil {
			return nil, err
		}
	}

	/* Read in the Frame Header. */
	frmHeaderWSQ, err := d.getCFrameHeaderWSQ()
	if err != nil {
		return nil, err
	}
	width := frmHeaderWSQ.width
	height := frmHeaderWSQ.height

	getCPpiWSQ()

	/* Build WSQ decomposition trees. */
	d.buildWSQTrees(width, height)

	/* Decode the Huffman encoded buffer blocks. */

	qdata, err := huffmanDecodeDataMem(&d, width*height)
	if err != nil {
		return nil, err
	}

	/* Decode the quantize wavelet subband buffer. */
	fdata, err := unquantize(&d, qdata, width, height)
	if err != nil {
		return nil, err
	}

	/* Done with quantized wavelet subband buffer. */
	//noinspection UnusedAssignment
	qdata = nil

	wsqReconstruct(&d, fdata, width, height)

	/* Convert floating point pixels to unsigned char pixels. */
	img := image.NewGray(image.Rect(0, 0, width, height))
	var idx int
	for r := 0; r < height; r++ {
		for c := 0; c < width; c++ {
			pixel := (fdata[idx] * frmHeaderWSQ.rScale) + frmHeaderWSQ.mShift
			pixel += 0.5

			if pixel < 0.0 {
				img.SetGray(c, r, color.Gray{Y: 0})
			} else if pixel > 255.0 {
				img.SetGray(c, r, color.Gray{Y: 255})
			} else {
				img.SetGray(c, r, color.Gray{Y: uint8(pixel)})
			}

			idx++
		}
	}

	return img, nil
}

func getCPpiWSQ() int {
	return -1
}

func huffmanDecodeDataMem(d *decoder, size int) ([]int, error) {
	qdata := make([]int, size)

	maxcode := make([]int, max_huffbits+1)
	mincode := make([]int, max_huffbits+1)
	valptr := make([]int, max_huffbits+1)

	ref, err := d.getCMarkerWSQ(tbls_n_sob)
	if err != nil {
		return nil, err
	}
	marker := &reference[int]{ref}

	bitCount := &reference[int]{0} /* bit count for getc_nextbits_wsq routine */
	nextByte := &reference[int]{0} /*next byte of buffer*/
	var hufftableId int            /* huffman table number */
	var ip int

	for marker.value != eoi_wsq {
		if marker.value != 0 {
			for marker.value != sob_wsq {
				if err := d.getCTableWSQ(marker.value); err != nil {
					return nil, err
				}
				marker.value, err = d.getCMarkerWSQ(tbls_n_sob)
				if err != nil {
					return nil, err
				}
				if marker.value == eoi_wsq {
					break
				}
			}
			if marker.value == eoi_wsq {
				break
			}
			hufftableId, err = d.getCBlockHeader() /* huffman table number */
			if err != nil {
				return nil, err
			}

			if d.tableDHT[hufftableId].tabdef != 1 {
				return nil, fmt.Errorf("ERROR : huffmanDecodeDataMem : huffman table undefined.")
			}

			/* the next two routines reconstruct the huffman tables */
			hufftable, _ := buildHuffsizes(d.tableDHT[hufftableId].huffbits, max_huffcounts_wsq)
			buildHuffcodes(hufftable)

			/* this routine builds a set of three tables used in decoding */
			/* the compressed buffer*/
			genDecodeTable(hufftable, maxcode, mincode, valptr, d.tableDHT[hufftableId].huffbits)

			bitCount.value = 0
			marker.value = 0
		}

		/* get next huffman category code from compressed input buffer stream */
		nodeptr, err := decodeDataMem(d, mincode, maxcode, valptr, d.tableDHT[hufftableId].huffvalues, bitCount, marker, nextByte)
		if err != nil {
			return nil, err
		}
		/* nodeptr  pointers for decoding */

		if nodeptr == -1 {
			continue
		}

		if nodeptr > 0 && nodeptr <= 100 {
			for n := 0; n < nodeptr; n++ {
				qdata[ip] = 0 /* z run */
				ip++
			}
		} else if nodeptr > 106 && nodeptr < 0xff {
			qdata[ip] = nodeptr - 180
			ip++
		} else if nodeptr == 101 {
			v, err := d.getCNextbitsWSQ(marker, bitCount, 8, nextByte)
			if err != nil {
				return nil, err
			}
			qdata[ip] = v
			ip++
		} else if nodeptr == 102 {
			v, err := d.getCNextbitsWSQ(marker, bitCount, 8, nextByte)
			if err != nil {
				return nil, err
			}
			qdata[ip] = -v
			ip++
		} else if nodeptr == 103 {
			v, err := d.getCNextbitsWSQ(marker, bitCount, 16, nextByte)
			if err != nil {
				return nil, err
			}
			qdata[ip] = v
			ip++
		} else if nodeptr == 104 {
			v, err := d.getCNextbitsWSQ(marker, bitCount, 16, nextByte)
			if err != nil {
				return nil, err
			}
			qdata[ip] = -v
			ip++
		} else if nodeptr == 105 {
			n, err := d.getCNextbitsWSQ(marker, bitCount, 8, nextByte)
			if err != nil {
				return nil, err
			}
			for n > 0 {
				n--
				qdata[ip] = 0
				ip++
			}
		} else if nodeptr == 106 {
			n, err := d.getCNextbitsWSQ(marker, bitCount, 16, nextByte)
			if err != nil {
				return nil, err
			}
			for n > 0 {
				n--
				qdata[ip] = 0
				ip++
			}
		} else {
			return nil, fmt.Errorf("ERROR: huffman_decode_data_mem : Invalid code (%d)", nodeptr)
		}
	}

	return qdata, nil
}

func buildHuffsizes(huffbits []int, maxHuffcounts int) ([]huffCode, int) {
	var huffcodeTable []huffCode /*table of huffman codes and sizes*/
	numberOfCodes := 1           /*the number codes for a given code size*/

	huffcodeTable = make([]huffCode, maxHuffcounts+1)

	var tempSize int

	for codeSize := 1; codeSize <= max_huffbits; codeSize++ {
		for numberOfCodes <= huffbits[codeSize-1] {
			huffcodeTable[tempSize] = huffCode{}
			huffcodeTable[tempSize].size = codeSize
			tempSize++
			numberOfCodes++
		}
		numberOfCodes = 1
	}

	huffcodeTable[tempSize] = huffCode{}
	huffcodeTable[tempSize].size = 0

	return huffcodeTable, tempSize
}

func buildHuffcodes(huffcodeTable []huffCode) {
	var tempCode int /*used to construct code word*/
	var pointer int  /*pointer to code word information*/

	tempSize := huffcodeTable[0].size
	if huffcodeTable[pointer].size == 0 {
		return
	}

	for {
		for {
			huffcodeTable[pointer].code = tempCode
			tempCode++
			pointer++
			if huffcodeTable[pointer].size != tempSize {
				break
			}
		}

		if huffcodeTable[pointer].size == 0 {
			return
		}

		for {
			tempCode <<= 1
			tempSize++
			if huffcodeTable[pointer].size == tempSize {
				break
			}
		}
		if huffcodeTable[pointer].size != tempSize {
			break
		}
	}
}

func genDecodeTable(huffcodeTable []huffCode, maxcode, mincode, valptr, huffbits []int) {
	for i := 0; i <= max_huffbits; i++ {
		maxcode[i] = 0
		mincode[i] = 0
		valptr[i] = 0
	}

	var i2 int
	for i := 1; i <= max_huffbits; i++ {
		if huffbits[i-1] == 0 {
			maxcode[i] = -1
			continue
		}
		valptr[i] = i2
		mincode[i] = huffcodeTable[i2].code
		i2 = i2 + huffbits[i-1] - 1
		maxcode[i] = huffcodeTable[i2].code
		i2++
	}
}

func decodeDataMem(d *decoder, mincode, maxcode, valptr, huffvalues []int, bitCount, marker, nextByte *reference[int]) (int, error) {
	code, err := d.getCNextbitsWSQ(marker, bitCount, 1, nextByte) /* becomes a huffman code word  (one bit at a time)*/
	if err != nil {
		return 0, err
	}
	if marker.value != 0 {
		return -1, nil
	}

	var inx int
	for inx = 1; code > maxcode[inx]; inx++ {
		tbits, err := d.getCNextbitsWSQ(marker, bitCount, 1, nextByte) /* becomes a huffman code word  (one bit at a time)*/
		if err != nil {
			return 0, err
		}
		code = ((code << 1) + tbits)

		if marker.value != 0 {
			return -1, nil
		}
	}

	inx2 := valptr[inx] + code - mincode[inx] /*increment variables*/
	return huffvalues[inx2], nil
}

func unquantize(d *decoder, sip []int, width, height int) ([]float32, error) {
	fip := make([]float32, width*height) /* floating point image */

	if d.tableDQT.dqtDef != 1 {
		return nil, fmt.Errorf("ERROR: unquantize : quantization table parameters not defined!")
	}

	binCenter := d.tableDQT.binCenter /* quantizer bin center */

	var sptr int
	for cnt := 0; cnt < num_subbands; cnt++ {
		if d.tableDQT.qBin[cnt] == 0.0 {
			continue
		}

		fptr := (d.qtree[cnt].y * width) + d.qtree[cnt].x

		for row := 0; row < d.qtree[cnt].leny; row++ {
			for col := 0; col < d.qtree[cnt].lenx; col++ {
				if sip[sptr] == 0 {
					fip[fptr] = 0
				} else if sip[sptr] > 0 {
					fip[fptr] = (d.tableDQT.qBin[cnt] * (float32(sip[sptr]) - binCenter)) + (d.tableDQT.zBin[cnt] / 2.0)
				} else if sip[sptr] < 0 {
					fip[fptr] = (d.tableDQT.qBin[cnt] * (float32(sip[sptr]) + binCenter)) - (d.tableDQT.zBin[cnt] / 2.0)
				} else {
					return nil, fmt.Errorf("ERROR : unquantize : invalid quantization pixel value")
				}
				fptr++
				sptr++
			}
			fptr += width - d.qtree[cnt].lenx
		}
	}

	return fip, nil
}

func wsqReconstruct(d *decoder, fdata []float32, width, height int) error {
	if d.tableDTT.lodef != 1 {
		return fmt.Errorf("ERROR: wsq_reconstruct : Lopass filter coefficients not defined")
	}

	if d.tableDTT.hidef != 1 {
		return fmt.Errorf("ERROR: wsq_reconstruct : Hipass filter coefficients not defined")
	}

	numPix := width * height
	/* Allocate temporary floating point pixmap. */
	fdataTemp := make([]float32, numPix)

	/* Reconstruct floating point pixmap from wavelet subband buffer. */
	for node := w_treelen - 1; node >= 0; node-- {
		fdataBse := (d.wtree[node].y * width) + d.wtree[node].x
		joinLets(fdataTemp, fdata, 0, fdataBse, d.wtree[node].lenx, d.wtree[node].leny,
			1, width,
			d.tableDTT.hifilt, d.tableDTT.hisz,
			d.tableDTT.lofilt, d.tableDTT.losz,
			d.wtree[node].invcl)
		joinLets(fdata, fdataTemp, fdataBse, 0, d.wtree[node].leny, d.wtree[node].lenx,
			width, 1,
			d.tableDTT.hifilt, d.tableDTT.hisz,
			d.tableDTT.lofilt, d.tableDTT.losz,
			d.wtree[node].invrw)
	}

	return nil
}

func joinLets(
	newdata,
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
	inv int) /* spectral inversion? */ {
	var lp0, lp1 int
	var hp0, hp1 int
	var lopass, hipass int /* lo/hi pass image pointers */
	var limg, himg int
	var pix, cl_rw int /* pixel counter and column/row counter */
	var i, da_ev int   /* if "scanline" is even or odd and */
	var loc, hoc int
	var hlen, llen int
	var nstr, pstr int
	var tap int
	var fi_ev int
	var olle, ohle, olre, ohre int
	var lle, lle2, lre, lre2 int
	var hle, hle2, hre, hre2 int
	var lpx, lspx int
	var lpxstr, lspxstr int
	var lstap, lotap int
	var hpx, hspx int
	var hpxstr, hspxstr int
	var hstap, hotap int
	var asym, fhre, ofhre int
	var ssfac, osfac, sfac float32

	da_ev = len2 % 2
	fi_ev = lsz % 2
	pstr = stride
	nstr = -pstr
	if da_ev != 0 {
		llen = (len2 + 1) / 2
		hlen = llen - 1
	} else {
		llen = len2 / 2
		hlen = llen
	}

	if fi_ev != 0 {
		asym = 0
		ssfac = 1.0
		ofhre = 0
		loc = (lsz - 1) / 4
		hoc = (hsz+1)/4 - 1
		lotap = ((lsz - 1) / 2) % 2
		hotap = ((hsz + 1) / 2) % 2
		if da_ev != 0 {
			olle = 0
			olre = 0
			ohle = 1
			ohre = 1
		} else {
			olle = 0
			olre = 1
			ohle = 1
			ohre = 0
		}
	} else {
		asym = 1
		ssfac = -1.0
		ofhre = 2
		loc = lsz/4 - 1
		hoc = hsz/4 - 1
		lotap = (lsz / 2) % 2
		hotap = (hsz / 2) % 2
		if da_ev != 0 {
			olle = 1
			olre = 0
			ohle = 1
			ohre = 1
		} else {
			olle = 1
			olre = 1
			ohle = 1
			ohre = 1
		}

		if loc == -1 {
			loc = 0
			olle = 0
		}
		if hoc == -1 {
			hoc = 0
			ohle = 0
		}

		for i = 0; i < hsz; i++ {
			hi[i] *= -1.0
		}
	}

	for cl_rw = 0; cl_rw < len1; cl_rw++ {
		limg = newIndex + cl_rw*pitch
		himg = limg
		newdata[himg] = 0.0
		newdata[himg+stride] = 0.0
		if inv != 0 {
			hipass = oldIndex + cl_rw*pitch
			lopass = hipass + stride*hlen
		} else {
			lopass = oldIndex + cl_rw*pitch
			hipass = lopass + stride*llen
		}

		lp0 = lopass
		lp1 = lp0 + (llen-1)*stride
		lspx = lp0 + (loc * stride)
		lspxstr = nstr
		lstap = lotap
		lle2 = olle
		lre2 = olre

		hp0 = hipass
		hp1 = hp0 + (hlen-1)*stride
		hspx = hp0 + (hoc * stride)
		hspxstr = nstr
		hstap = hotap
		hle2 = ohle
		hre2 = ohre
		osfac = ssfac

		for pix = 0; pix < hlen; pix++ {
			for tap = lstap; tap >= 0; tap-- {
				lle = lle2
				lre = lre2
				lpx = lspx
				lpxstr = lspxstr

				newdata[limg] = olddata[lpx] * lo[tap]
				for i = tap + 2; i < lsz; i += 2 {
					if lpx == lp0 {
						if lle != 0 {
							lpxstr = 0
							lle = 0
						} else {
							lpxstr = pstr
						}
					}
					if lpx == lp1 {
						if lre != 0 {
							lpxstr = 0
							lre = 0
						} else {
							lpxstr = nstr
						}
					}
					lpx += lpxstr
					newdata[limg] += olddata[lpx] * lo[i]
				}
				limg += stride
			}
			if lspx == lp0 {
				if lle2 != 0 {
					lspxstr = 0
					lle2 = 0
				} else {
					lspxstr = pstr
				}
			}
			lspx += lspxstr
			lstap = 1

			for tap = hstap; tap >= 0; tap-- {
				hle = hle2
				hre = hre2
				hpx = hspx
				hpxstr = hspxstr
				fhre = ofhre
				sfac = osfac

				for i = tap; i < hsz; i += 2 {
					if hpx == hp0 {
						if hle != 0 {
							hpxstr = 0
							hle = 0
						} else {
							hpxstr = pstr
							sfac = 1.0
						}
					}
					if hpx == hp1 {
						if hre != 0 {
							hpxstr = 0
							hre = 0
							if asym != 0 && da_ev != 0 {
								hre = 1
								fhre--
								sfac = float32(fhre)
								if sfac == 0.0 {
									hre = 0
								}
							}
						} else {
							hpxstr = nstr
							if asym != 0 {
								sfac = -1.0
							}
						}
					}
					newdata[himg] += olddata[hpx] * hi[i] * sfac
					hpx += hpxstr
				}
				himg += stride
			}
			if hspx == hp0 {
				if hle2 != 0 {
					hspxstr = 0
					hle2 = 0
				} else {
					hspxstr = pstr
					osfac = 1.0
				}
			}
			hspx += hspxstr
			hstap = 1
		}

		if da_ev != 0 {
			if lotap != 0 {
				lstap = 1
			} else {
				lstap = 0
			}
		} else if lotap != 0 {
			lstap = 2
		} else {
			lstap = 1
		}

		for tap = 1; tap >= lstap; tap-- {
			lle = lle2
			lre = lre2
			lpx = lspx
			lpxstr = lspxstr

			newdata[limg] = olddata[lpx] * lo[tap]
			for i = tap + 2; i < lsz; i += 2 {
				if lpx == lp0 {
					if lle != 0 {
						lpxstr = 0
						lle = 0
					} else {
						lpxstr = pstr
					}
				}
				if lpx == lp1 {
					if lre != 0 {
						lpxstr = 0
						lre = 0
					} else {
						lpxstr = nstr
					}
				}
				lpx += lpxstr
				newdata[limg] += olddata[lpx] * lo[i]
			}
			limg += stride
		}

		if da_ev != 0 {
			if hotap != 0 {
				hstap = 1
			} else {
				hstap = 0
			}
			if hsz == 2 {
				hspx -= hspxstr
				fhre = 1
			}
		} else if hotap != 0 {
			hstap = 2
		} else {
			hstap = 1
		}

		for tap = 1; tap >= hstap; tap-- {
			hle = hle2
			hre = hre2
			hpx = hspx
			hpxstr = hspxstr
			sfac = osfac
			if hsz != 2 {
				fhre = ofhre
			}

			for i = tap; i < hsz; i += 2 {
				if hpx == hp0 {
					if hle != 0 {
						hpxstr = 0
						hle = 0
					} else {
						hpxstr = pstr
						sfac = 1.0
					}
				}
				if hpx == hp1 {
					if hre != 0 {
						hpxstr = 0
						hre = 0
						if asym != 0 && da_ev != 0 {
							hre = 1
							fhre--
							sfac = float32(fhre)
							if sfac == 0.0 {
								hre = 0
							}
						}
					} else {
						hpxstr = nstr
						if asym != 0 {
							sfac = -1.0
						}
					}
				}
				newdata[himg] += olddata[hpx] * hi[i] * sfac
				hpx += hpxstr
			}
			himg += stride
		}
	}

	if fi_ev == 0 {
		for i = 0; i < hsz; i++ {
			hi[i] *= -1.0
		}
	}
}
