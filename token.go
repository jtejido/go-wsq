package wsq

import (
	"io"
)

const (
	max_dht_tables     = 8
	max_subbands       = 64
	max_huffbits       = 16
	max_huffcounts_wsq = 256
	max_huffcoeff      = 74 /* -73 .. +74 */
	max_huffzrun       = 100

	max_hifilt = 7
	max_lofilt = 9
	w_treelen  = 20
	q_treelen  = 64
	int32size  = 4
	int16size  = 2
	/* wsq marker definitions */
	soi_wsq = 0xffa0
	eoi_wsq = 0xffa1
	sof_wsq = 0xffa2
	sob_wsq = 0xffa3
	dtt_wsq = 0xffa4
	dqt_wsq = 0xffa5
	dht_wsq = 0xffa6
	drt_wsq = 0xffa7
	com_wsq = 0xffa8

	strt_subband_2     = 19
	strt_subband_3     = 52
	num_subbands       = 60
	strt_subband_del   = num_subbands
	strt_size_region_2 = 4
	strt_size_region_3 = 51

	coeff_code = 0
	run_code   = 1

	variance_thresh float32 = 1.01

	/* case for getting any marker. */
	any_wsq    = 0xffff
	tbls_n_sof = 2
	tbls_n_sob = tbls_n_sof + 2
)

var (
	bitmask            = []int{0x00, 0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff}
	hi_filt_even_8x8_1 = []float32{
		0.03226944131446922,
		-0.05261415011924844,
		-0.18870142780632693,
		0.60328894481393847,
		-0.60328894481393847,
		0.18870142780632693,
		0.05261415011924844,
		-0.03226944131446922,
	}

	lo_filt_even_8x8_1 = []float32{
		0.07565691101399093,
		-0.12335584105275092,
		-0.09789296778409587,
		0.85269867900940344,
		0.85269867900940344,
		-0.09789296778409587,
		-0.12335584105275092,
		0.07565691101399093,
	}

	hi_filt_not_even_8x8_1 = []float32{
		0.06453888262893845,
		-0.04068941760955844,
		-0.41809227322221221,
		0.78848561640566439,
		-0.41809227322221221,
		-0.04068941760955844,
		0.06453888262893845,
	}

	lo_filt_not_even_8x8_1 = []float32{
		0.03782845550699546,
		-0.02384946501938000,
		-0.11062440441842342,
		0.37740285561265380,
		0.85269867900940344,
		0.37740285561265380,
		-0.11062440441842342,
		-0.02384946501938000,
		0.03782845550699546,
	}
)

type tableDTT struct {
	lofilt                   []float32
	hifilt                   []float32
	losz, hisz, lodef, hidef int
}
type tableDQT struct {
	binCenter float32
	qBin      [max_subbands]float32
	zBin      [max_subbands]float32
	dqtDef    rune
}

type huffmanTable struct {
	tableLen             int
	bytesLeft            int
	tableId              int
	huffbits, huffvalues []int
}

type tableDHT struct {
	tabdef     byte
	huffbits   []int // MAX_HUFFBITS
	huffvalues []int // MAX_HUFFCOUNTS_WSQ + 1
}
type wavletTree struct {
	x, y, lenx, leny, invrw, invcl int
}

type quantTree struct {
	x, y, lenx, leny int
}

type quantization struct {
	q        float32 /* quantization level */
	cr       float32 /* compression ratio */
	r        float32 /* compression bitrate */
	qbss_t   [max_subbands]float32
	qbss     [max_subbands]float32
	qzbs     [max_subbands]float32
	variance [max_subbands]float32
}

type ref interface {
	int | float32 | bool | []int
}

type reference[V ref] struct {
	value V
}

type headerFrm struct {
	black, white, width, height int
	mShift, rScale              float32
	wsqEncoder, software        int
}

type huffCode struct {
	size int
	code int
}

func readByte(r io.ByteReader) (byte, error) {
	b, err := r.ReadByte()
	if err == io.EOF {
		err = io.ErrUnexpectedEOF
	}
	return b, err
}

type token struct {
	tableDHT []*tableDHT
	tableDTT tableDTT
	tableDQT tableDQT
	wtree    []wavletTree
	qtree    []quantTree
}

func (t *token) buildWSQTrees(width, height int) {
	/* Build a W-TREE structure for the image. */
	t.buildWTree(w_treelen, width, height)
	/* Build a Q-TREE structure for the image. */
	t.buildQTree(q_treelen)
}

func (t *token) buildWTree(wtreelen, width, height int) {
	var lenx, lenx2, leny, leny2 int /* starting lengths of sections of the image being split into subbands */
	t.wtree = make([]wavletTree, wtreelen)
	for i := 0; i < wtreelen; i++ {
		t.wtree[i] = wavletTree{
			invrw: 0,
			invcl: 0,
		}
	}

	t.wtree[2].invrw = 1
	t.wtree[4].invrw = 1
	t.wtree[7].invrw = 1
	t.wtree[9].invrw = 1
	t.wtree[11].invrw = 1
	t.wtree[13].invrw = 1
	t.wtree[16].invrw = 1
	t.wtree[18].invrw = 1
	t.wtree[3].invcl = 1
	t.wtree[5].invcl = 1
	t.wtree[8].invcl = 1
	t.wtree[9].invcl = 1
	t.wtree[12].invcl = 1
	t.wtree[13].invcl = 1
	t.wtree[17].invcl = 1
	t.wtree[18].invcl = 1

	t.wtree4(0, 1, width, height, 0, 0, 1)

	if (t.wtree[1].lenx % 2) == 0 {
		lenx = t.wtree[1].lenx / 2
		lenx2 = lenx
	} else {
		lenx = (t.wtree[1].lenx + 1) / 2
		lenx2 = lenx - 1
	}

	if (t.wtree[1].leny % 2) == 0 {
		leny = t.wtree[1].leny / 2
		leny2 = leny
	} else {
		leny = (t.wtree[1].leny + 1) / 2
		leny2 = leny - 1
	}

	t.wtree4(4, 6, lenx2, leny, lenx, 0, 0)
	t.wtree4(5, 10, lenx, leny2, 0, leny, 0)
	t.wtree4(14, 15, lenx, leny, 0, 0, 0)

	t.wtree[19].x = 0
	t.wtree[19].y = 0
	if (t.wtree[15].lenx % 2) == 0 {
		t.wtree[19].lenx = t.wtree[15].lenx / 2
	} else {
		t.wtree[19].lenx = (t.wtree[15].lenx + 1) / 2
	}

	if (t.wtree[15].leny % 2) == 0 {
		t.wtree[19].leny = t.wtree[15].leny / 2
	} else {
		t.wtree[19].leny = (t.wtree[15].leny + 1) / 2
	}
}

func (t *token) wtree4(start1, start2, lenx, leny, x, y, stop1 int) {
	var evenx, eveny int /* Check length of subband for even or odd */
	var p1, p2 int       /* w_tree locations for storing subband sizes and locations */

	p1 = start1
	p2 = start2

	evenx = lenx % 2
	eveny = leny % 2

	t.wtree[p1].x = x
	t.wtree[p1].y = y
	t.wtree[p1].lenx = lenx
	t.wtree[p1].leny = leny

	t.wtree[p2].x = x
	t.wtree[p2+2].x = x
	t.wtree[p2].y = y
	t.wtree[p2+1].y = y

	if evenx == 0 {
		t.wtree[p2].lenx = lenx / 2
		t.wtree[p2+1].lenx = t.wtree[p2].lenx
	} else {
		if p1 == 4 {
			t.wtree[p2].lenx = (lenx - 1) / 2
			t.wtree[p2+1].lenx = t.wtree[p2].lenx + 1
		} else {
			t.wtree[p2].lenx = (lenx + 1) / 2
			t.wtree[p2+1].lenx = t.wtree[p2].lenx - 1
		}
	}
	t.wtree[p2+1].x = t.wtree[p2].lenx + x
	if stop1 == 0 {
		t.wtree[p2+3].lenx = t.wtree[p2+1].lenx
		t.wtree[p2+3].x = t.wtree[p2+1].x
	}
	t.wtree[p2+2].lenx = t.wtree[p2].lenx

	if eveny == 0 {
		t.wtree[p2].leny = leny / 2
		t.wtree[p2+2].leny = t.wtree[p2].leny
	} else {
		if p1 == 5 {
			t.wtree[p2].leny = (leny - 1) / 2
			t.wtree[p2+2].leny = t.wtree[p2].leny + 1
		} else {
			t.wtree[p2].leny = (leny + 1) / 2
			t.wtree[p2+2].leny = t.wtree[p2].leny - 1
		}
	}
	t.wtree[p2+2].y = t.wtree[p2].leny + y
	if stop1 == 0 {
		t.wtree[p2+3].leny = t.wtree[p2+2].leny
		t.wtree[p2+3].y = t.wtree[p2+2].y
	}
	t.wtree[p2+1].leny = t.wtree[p2].leny
}

func (t *token) buildQTree(qtreelen int) {
	t.qtree = make([]quantTree, qtreelen)
	for i := 0; i < len(t.qtree); i++ {
		t.qtree[i] = quantTree{}
	}

	t.qtree16(3, t.wtree[14].lenx, t.wtree[14].leny, t.wtree[14].x, t.wtree[14].y, 0, 0)
	t.qtree16(19, t.wtree[4].lenx, t.wtree[4].leny, t.wtree[4].x, t.wtree[4].y, 0, 1)
	t.qtree16(48, t.wtree[0].lenx, t.wtree[0].leny, t.wtree[0].x, t.wtree[0].y, 0, 0)
	t.qtree16(35, t.wtree[5].lenx, t.wtree[5].leny, t.wtree[5].x, t.wtree[5].y, 1, 0)
	t.qtree4(0, t.wtree[19].lenx, t.wtree[19].leny, t.wtree[19].x, t.wtree[19].y)
}

func (t *token) qtree16(start, lenx, leny, x, y, rw, cl int) {
	var tempx, temp2x int /* temporary x values */
	var tempy, temp2y int /* temporary y values */
	var evenx, eveny int  /* Check length of subband for even or odd */
	var p int             /* indicates subband information being stored */

	p = start
	evenx = lenx % 2
	eveny = leny % 2

	if evenx == 0 {
		tempx = lenx / 2
		temp2x = tempx
	} else {
		if cl != 0 {
			temp2x = (lenx + 1) / 2
			tempx = temp2x - 1
		} else {
			tempx = (lenx + 1) / 2
			temp2x = tempx - 1
		}
	}

	if eveny == 0 {
		tempy = leny / 2
		temp2y = tempy
	} else {
		if rw != 0 {
			temp2y = (leny + 1) / 2
			tempy = temp2y - 1
		} else {
			tempy = (leny + 1) / 2
			temp2y = tempy - 1
		}
	}

	evenx = tempx % 2
	eveny = tempy % 2

	t.qtree[p].x = x
	t.qtree[p+2].x = x
	t.qtree[p].y = y
	t.qtree[p+1].y = y
	if evenx == 0 {
		t.qtree[p].lenx = tempx / 2
		t.qtree[p+1].lenx = t.qtree[p].lenx
		t.qtree[p+2].lenx = t.qtree[p].lenx
		t.qtree[p+3].lenx = t.qtree[p].lenx
	} else {
		t.qtree[p].lenx = (tempx + 1) / 2
		t.qtree[p+1].lenx = t.qtree[p].lenx - 1
		t.qtree[p+2].lenx = t.qtree[p].lenx
		t.qtree[p+3].lenx = t.qtree[p+1].lenx
	}
	t.qtree[p+1].x = x + t.qtree[p].lenx
	t.qtree[p+3].x = t.qtree[p+1].x
	if eveny == 0 {
		t.qtree[p].leny = tempy / 2
		t.qtree[p+1].leny = t.qtree[p].leny
		t.qtree[p+2].leny = t.qtree[p].leny
		t.qtree[p+3].leny = t.qtree[p].leny
	} else {
		t.qtree[p].leny = (tempy + 1) / 2
		t.qtree[p+1].leny = t.qtree[p].leny
		t.qtree[p+2].leny = t.qtree[p].leny - 1
		t.qtree[p+3].leny = t.qtree[p+2].leny
	}
	t.qtree[p+2].y = y + t.qtree[p].leny
	t.qtree[p+3].y = t.qtree[p+2].y

	evenx = temp2x % 2

	t.qtree[p+4].x = x + tempx
	t.qtree[p+6].x = t.qtree[p+4].x
	t.qtree[p+4].y = y
	t.qtree[p+5].y = y
	t.qtree[p+6].y = t.qtree[p+2].y
	t.qtree[p+7].y = t.qtree[p+2].y
	t.qtree[p+4].leny = t.qtree[p].leny
	t.qtree[p+5].leny = t.qtree[p].leny
	t.qtree[p+6].leny = t.qtree[p+2].leny
	t.qtree[p+7].leny = t.qtree[p+2].leny
	if evenx == 0 {
		t.qtree[p+4].lenx = temp2x / 2
		t.qtree[p+5].lenx = t.qtree[p+4].lenx
		t.qtree[p+6].lenx = t.qtree[p+4].lenx
		t.qtree[p+7].lenx = t.qtree[p+4].lenx
	} else {
		t.qtree[p+5].lenx = (temp2x + 1) / 2
		t.qtree[p+4].lenx = t.qtree[p+5].lenx - 1
		t.qtree[p+6].lenx = t.qtree[p+4].lenx
		t.qtree[p+7].lenx = t.qtree[p+5].lenx
	}
	t.qtree[p+5].x = t.qtree[p+4].x + t.qtree[p+4].lenx
	t.qtree[p+7].x = t.qtree[p+5].x

	eveny = temp2y % 2

	t.qtree[p+8].x = x
	t.qtree[p+9].x = t.qtree[p+1].x
	t.qtree[p+10].x = x
	t.qtree[p+11].x = t.qtree[p+1].x
	t.qtree[p+8].y = y + tempy
	t.qtree[p+9].y = t.qtree[p+8].y
	t.qtree[p+8].lenx = t.qtree[p].lenx
	t.qtree[p+9].lenx = t.qtree[p+1].lenx
	t.qtree[p+10].lenx = t.qtree[p].lenx
	t.qtree[p+11].lenx = t.qtree[p+1].lenx
	if eveny == 0 {
		t.qtree[p+8].leny = temp2y / 2
		t.qtree[p+9].leny = t.qtree[p+8].leny
		t.qtree[p+10].leny = t.qtree[p+8].leny
		t.qtree[p+11].leny = t.qtree[p+8].leny
	} else {
		t.qtree[p+10].leny = (temp2y + 1) / 2
		t.qtree[p+11].leny = t.qtree[p+10].leny
		t.qtree[p+8].leny = t.qtree[p+10].leny - 1
		t.qtree[p+9].leny = t.qtree[p+8].leny
	}
	t.qtree[p+10].y = t.qtree[p+8].y + t.qtree[p+8].leny
	t.qtree[p+11].y = t.qtree[p+10].y

	t.qtree[p+12].x = t.qtree[p+4].x
	t.qtree[p+13].x = t.qtree[p+5].x
	t.qtree[p+14].x = t.qtree[p+4].x
	t.qtree[p+15].x = t.qtree[p+5].x
	t.qtree[p+12].y = t.qtree[p+8].y
	t.qtree[p+13].y = t.qtree[p+8].y
	t.qtree[p+14].y = t.qtree[p+10].y
	t.qtree[p+15].y = t.qtree[p+10].y
	t.qtree[p+12].lenx = t.qtree[p+4].lenx
	t.qtree[p+13].lenx = t.qtree[p+5].lenx
	t.qtree[p+14].lenx = t.qtree[p+4].lenx
	t.qtree[p+15].lenx = t.qtree[p+5].lenx
	t.qtree[p+12].leny = t.qtree[p+8].leny
	t.qtree[p+13].leny = t.qtree[p+8].leny
	t.qtree[p+14].leny = t.qtree[p+10].leny
	t.qtree[p+15].leny = t.qtree[p+10].leny
}

func (t *token) qtree4(start, lenx, leny, x, y int) {
	var evenx, eveny int /* Check length of subband for even or odd */
	var p int            /* indicates subband information being stored */

	p = start
	evenx = lenx % 2
	eveny = leny % 2

	t.qtree[p].x = x
	t.qtree[p+2].x = x
	t.qtree[p].y = y
	t.qtree[p+1].y = y
	if evenx == 0 {
		t.qtree[p].lenx = lenx / 2
		t.qtree[p+1].lenx = t.qtree[p].lenx
		t.qtree[p+2].lenx = t.qtree[p].lenx
		t.qtree[p+3].lenx = t.qtree[p].lenx
	} else {
		t.qtree[p].lenx = (lenx + 1) / 2
		t.qtree[p+1].lenx = t.qtree[p].lenx - 1
		t.qtree[p+2].lenx = t.qtree[p].lenx
		t.qtree[p+3].lenx = t.qtree[p+1].lenx
	}
	t.qtree[p+1].x = x + t.qtree[p].lenx
	t.qtree[p+3].x = t.qtree[p+1].x
	if eveny == 0 {
		t.qtree[p].leny = leny / 2
		t.qtree[p+1].leny = t.qtree[p].leny
		t.qtree[p+2].leny = t.qtree[p].leny
		t.qtree[p+3].leny = t.qtree[p].leny
	} else {
		t.qtree[p].leny = (leny + 1) / 2
		t.qtree[p+1].leny = t.qtree[p].leny
		t.qtree[p+2].leny = t.qtree[p].leny - 1
		t.qtree[p+3].leny = t.qtree[p+2].leny
	}
	t.qtree[p+2].y = y + t.qtree[p].leny
	t.qtree[p+3].y = t.qtree[p+2].y
}
