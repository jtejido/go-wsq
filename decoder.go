package wsq

import (
	"fmt"
	"io"
)

type decoder struct {
	token
	r        io.Reader
	pointer  int
	comments string
}

func (d *decoder) init() {
	d.tableDHT = make([]*tableDHT, max_dht_tables)

	for i := 0; i < max_dht_tables; i++ {
		d.tableDHT[i] = new(tableDHT)
		d.tableDHT[i].tabdef = 0
		d.tableDHT[i].huffbits = make([]int, max_huffbits)
		d.tableDHT[i].huffvalues = make([]int, max_huffcounts_wsq+1)
	}
}

func (d *decoder) getCBlockHeader() (int, error) {
	_, err := d.readShort() /* block header size */
	if err != nil {
		return 0, err
	}
	return d.readByte()
}

func (d *decoder) getCTableWSQ(marker int) (err error) {
	switch marker {
	case dtt_wsq:
		d.getCTransformTable()
		return nil
	case dqt_wsq:
		d.getCQuantizationTable()
		return nil
	case dht_wsq:
		d.getCHuffmanTableWSQ()
		return nil
	case com_wsq:
		//shams: i don't use return value
		d.comments, err = d.getCComment()
		if err != nil {
			return err
		}
		return nil
	default:
		return fmt.Errorf("ERROR: getCTableWSQ : Invalid table defined : %d", marker)
	}
}

func (d *decoder) getCFrameHeaderWSQ() (headerFrm headerFrm, err error) {
	_, err = d.readShort() /* header size */
	if err != nil {
		return
	}

	headerFrm.black, err = d.readByte()
	if err != nil {
		return
	}

	headerFrm.white, err = d.readByte()
	if err != nil {
		return
	}

	headerFrm.height, err = d.readShort()
	if err != nil {
		return
	}

	headerFrm.width, err = d.readShort()
	if err != nil {
		return
	}
	var scale, shrtDat int
	scale, err = d.readByte() /* exponent scaling parameter */
	if err != nil {
		return
	}

	shrtDat, err = d.readShort() /* buffer pointer */
	if err != nil {
		return
	}

	headerFrm.mShift = float32(shrtDat)
	for scale > 0 {
		headerFrm.mShift /= 10.0
		scale--
	}

	scale, err = d.readByte()
	if err != nil {
		return
	}
	shrtDat, err = d.readShort()
	if err != nil {
		return
	}
	headerFrm.rScale = float32(shrtDat)
	for scale > 0 {
		headerFrm.rScale /= 10.0
		scale--
	}

	headerFrm.wsqEncoder, err = d.readByte()
	if err != nil {
		return
	}
	headerFrm.software, err = d.readShort()
	if err != nil {
		return
	}

	return
}

func (d *decoder) getCComment() (string, error) {
	s, err := d.readShort()
	if err != nil {
		return "", err
	}
	size := s - 2
	comment, err := d.readBytes(size)
	return string(comment), err
}

func (d *decoder) getCTransformTable() (err error) {
	// read header Size;
	d.readShort()

	d.tableDTT.hisz, err = d.readByte()
	if err != nil {
		return
	}
	d.tableDTT.losz, err = d.readByte()
	if err != nil {
		return
	}

	d.tableDTT.hifilt = make([]float32, d.tableDTT.hisz)
	d.tableDTT.lofilt = make([]float32, d.tableDTT.losz)

	var aSize int
	if d.tableDTT.hisz%2 != 0 {
		aSize = (d.tableDTT.hisz + 1) / 2
	} else {
		aSize = d.tableDTT.hisz / 2
	}

	aLofilt := make([]float32, aSize)

	aSize--
	for cnt := 0; cnt <= aSize; cnt++ {
		sign, err := d.readByte()
		if err != nil {
			return err
		}
		scale, err := d.readByte()
		if err != nil {
			return err
		}
		shrtDat, err := d.readLong()
		if err != nil {
			return err
		}
		aLofilt[cnt] = float32(shrtDat)

		for scale > 0 {
			aLofilt[cnt] /= 10.0
			scale--
		}

		if sign != 0 {
			aLofilt[cnt] *= -1.0
		}

		if d.tableDTT.hisz%2 != 0 {
			d.tableDTT.hifilt[cnt+aSize] = float32(intSign(cnt)) * aLofilt[cnt]
			if cnt > 0 {
				d.tableDTT.hifilt[aSize-cnt] = d.tableDTT.hifilt[cnt+aSize]
			}
		} else {
			d.tableDTT.hifilt[cnt+aSize+1] = float32(intSign(cnt)) * aLofilt[cnt]
			d.tableDTT.hifilt[aSize-cnt] = -1 * d.tableDTT.hifilt[cnt+aSize+1]
		}
	}

	if d.tableDTT.losz%2 != 0 {
		aSize = (d.tableDTT.losz + 1) / 2
	} else {
		aSize = d.tableDTT.losz / 2
	}

	aHifilt := make([]float32, aSize)

	aSize--
	for cnt := 0; cnt <= aSize; cnt++ {
		sign, err := d.readByte()
		if err != nil {
			return err
		}
		scale, err := d.readByte()
		if err != nil {
			return err
		}
		shrtDat, err := d.readLong()
		if err != nil {
			return err
		}

		aHifilt[cnt] = float32(shrtDat)

		for scale > 0 {
			aHifilt[cnt] /= 10.0
			scale--
		}

		if sign != 0 {
			aHifilt[cnt] *= -1.0
		}

		if d.tableDTT.losz%2 != 0 {
			d.tableDTT.lofilt[cnt+aSize] = float32(intSign(cnt)) * aHifilt[cnt]
			if cnt > 0 {
				d.tableDTT.lofilt[aSize-cnt] = d.tableDTT.lofilt[cnt+aSize]
			}
		} else {
			d.tableDTT.lofilt[cnt+aSize+1] = float32(intSign(cnt+1)) * aHifilt[cnt]
			d.tableDTT.lofilt[aSize-cnt] = d.tableDTT.lofilt[cnt+aSize+1]
		}
	}

	d.tableDTT.lodef = 1
	d.tableDTT.hidef = 1
	return nil
}

func (d *decoder) getCQuantizationTable() (err error) {
	d.readShort()              /* header size */
	scale, err := d.readByte() /* scaling parameter */
	if err != nil {
		return err
	}
	shrtDat, err := d.readShort() /* counter and temp short buffer */
	if err != nil {
		return err
	}

	d.tableDQT.binCenter = float32(shrtDat)
	for scale > 0 {
		d.tableDQT.binCenter /= 10.0
		scale--
	}

	for cnt := 0; cnt < max_subbands; cnt++ {
		scale, err = d.readByte()
		if err != nil {
			return err
		}
		shrtDat, err = d.readShort()
		if err != nil {
			return err
		}
		d.tableDQT.qBin[cnt] = float32(shrtDat)
		for scale > 0 {
			d.tableDQT.qBin[cnt] /= 10.0
			scale--
		}

		scale, err = d.readByte()
		if err != nil {
			return err
		}
		shrtDat, err = d.readShort()
		if err != nil {
			return err
		}
		d.tableDQT.zBin[cnt] = float32(shrtDat)
		for scale > 0 {
			d.tableDQT.zBin[cnt] /= 10.0
			scale--
		}
	}

	d.tableDQT.dqtDef = 1
	return nil
}

func (d *decoder) getCHuffmanTableWSQ() error {
	/* First time, read table len. */
	firstHuffmanTable, err := d.getCHuffmanTable(max_huffcounts_wsq, 0, true)
	if err != nil {
		return err
	}

	/* Store table into global structure list. */
	tableId := firstHuffmanTable.tableId
	copy(d.tableDHT[tableId].huffbits, firstHuffmanTable.huffbits)
	copy(d.tableDHT[tableId].huffvalues, firstHuffmanTable.huffvalues)
	d.tableDHT[tableId].tabdef = 1

	bytesLeft := firstHuffmanTable.bytesLeft
	for bytesLeft != 0 {
		/* Read next table without reading table len. */
		huffmantable, err := d.getCHuffmanTable(max_huffcounts_wsq, bytesLeft, false)
		if err != nil {
			return err
		}

		/* If table is already defined ... */
		tableId = huffmantable.tableId
		if d.tableDHT[tableId].tabdef != 0 {
			return fmt.Errorf("ERROR : getCHuffmanTableWSQ : huffman table already defined.")
		}

		/* Store table into global structure list. */
		copy(d.tableDHT[tableId].huffbits, huffmantable.huffbits)
		copy(d.tableDHT[tableId].huffvalues, huffmantable.huffvalues)
		d.tableDHT[tableId].tabdef = 1
		bytesLeft = huffmantable.bytesLeft
	}

	return nil
}

func (d *decoder) getCHuffmanTable(maxHuffcounts, bytesLeft int, readTableLen bool) (htable *huffmanTable, err error) {
	htable = &huffmanTable{}

	/* table_len */
	if readTableLen {
		htable.tableLen, err = d.readShort()
		if err != nil {
			return nil, err
		}
		htable.bytesLeft = htable.tableLen - 2
		bytesLeft = htable.bytesLeft
	} else {
		htable.bytesLeft = bytesLeft
	}

	/* If no bytes left ... */
	if bytesLeft <= 0 {
		return nil, fmt.Errorf("ERROR : getCHuffmanTable : no huffman table bytes remaining")
	}

	/* Table ID */
	htable.tableId, err = d.readByte()
	if err != nil {
		return nil, err
	}
	htable.bytesLeft--

	htable.huffbits = make([]int, max_huffbits)
	var numHufvals int
	/* L1 ... L16 */
	for i := 0; i < max_huffbits; i++ {
		htable.huffbits[i], err = d.readByte()
		if err != nil {
			return nil, err
		}
		numHufvals += htable.huffbits[i]
	}
	htable.bytesLeft -= max_huffbits

	if numHufvals > maxHuffcounts+1 {
		return nil, fmt.Errorf("ERROR : getCHuffmanTable : numHufvals is larger than MAX_HUFFCOUNTS")
	}

	/* Could allocate only the amount needed ... then we wouldn't */
	/* need to pass MAX_HUFFCOUNTS. */
	htable.huffvalues = make([]int, maxHuffcounts+1)

	/* V1,1 ... V16,16 */
	for i := 0; i < numHufvals; i++ {
		htable.huffvalues[i], err = d.readByte()
		if err != nil {
			return nil, err
		}
	}
	htable.bytesLeft -= numHufvals

	return
}

func (d *decoder) getCMarkerWSQ(tt int) (int, error) {
	marker, err := d.readShort()
	if err != nil {
		return 0, err
	}

	switch tt {
	case soi_wsq:
		if marker != soi_wsq {
			return 0, fmt.Errorf("ERROR : getCMarkerWSQ : No SOI marker : %d", marker)
		}

		return marker, nil

	case tbls_n_sof:
		if marker != dtt_wsq && marker != dqt_wsq && marker != dht_wsq && marker != sof_wsq && marker != com_wsq && marker != eoi_wsq {
			return 0, fmt.Errorf("ERROR : getc_marker_wsq : No SOF, Table, or comment markers : %d", marker)
		}

		return marker, nil

	case tbls_n_sob:
		if marker != dtt_wsq && marker != dqt_wsq && marker != dht_wsq && marker != sob_wsq && marker != com_wsq && marker != eoi_wsq {
			return 0, fmt.Errorf("ERROR : getc_marker_wsq : No SOB, Table, or comment markers :  %d", marker)
		}
		return marker, nil
	case any_wsq:
		if (marker & 0xff00) != 0xff00 {
			return 0, fmt.Errorf("ERROR : getc_marker_wsq : no marker found : %d", marker)
		}

		/* Added by MDG on 03-07-05 */
		if (marker < soi_wsq) || (marker > com_wsq) {
			return 0, fmt.Errorf("ERROR : getc_marker_wsq : not a valid marker : %d", marker)
		}

		return marker, nil
	default:
		return 0, fmt.Errorf("ERROR : getc_marker_wsq : Invalid marker : %d", marker)
	}
}

func (d *decoder) getCNextbitsWSQ(marker, bitCount *reference[int], bitsReq int, nextByte *reference[int]) (bits int, err error) {
	if bitCount.value == 0 {
		nextByte.value, err = d.readByte()
		if err != nil {
			return
		}

		bitCount.value = 8
		if nextByte.value == 0xFF {
			var code2 int
			code2, err = d.readByte() /*stuffed byte of buffer*/
			if err != nil {
				return
			}
			if code2 != 0x00 && bitsReq == 1 {
				marker.value = (nextByte.value << 8) | code2
				return 1, nil
			}
			if code2 != 0x00 {
				return 0, fmt.Errorf("ERROR: getCNextbitsWSQ : No stuffed zeros.")
			}
		}
	}

	var tbits int      /*bits of current buffer byte requested*/
	var bitsNeeded int /*additional bits required to finish request*/
	if bitsReq <= bitCount.value {
		bits = (nextByte.value >> (bitCount.value - bitsReq)) & (bitmask[bitsReq])
		bitCount.value -= bitsReq
		nextByte.value &= bitmask[bitCount.value]
	} else {
		bitsNeeded = bitsReq - bitCount.value /*additional bits required to finish request*/
		bits = nextByte.value << bitsNeeded
		bitCount.value = 0
		tbits, err = d.getCNextbitsWSQ(marker, bitCount, bitsNeeded, nextByte)
		if err != nil {
			return 0, err
		}
		bits |= tbits
	}

	return bits, nil
}

func (d *decoder) Read(p []byte) (n int, err error) {
	n, err = d.r.Read(p)
	d.pointer += n

	return n, err
}

// read wraps Read returning only error.
func (d *decoder) read(p []byte) (err error) {
	_, err = d.Read(p)
	return err
}

func allocateBytes(size int) []byte {
	return make([]byte, size)
}

func (d *decoder) readLong() (int64, error) {
	b := allocateBytes(int32size)
	if err := d.read(b); err != nil { // read required bytes amount counting taken bytes internally
		return 0, err
	}
	b1 := int64(b[0])
	b2 := int64(b[1])
	b3 := int64(b[2])
	b4 := int64(b[3])
	value := (0xff&b1)<<24 | (0xff&b2)<<16 | (0xff&b3)<<8 | (0xff & b4)
	return value, nil
}

func (d *decoder) readShort() (int, error) {
	b := allocateBytes(int16size)
	if err := d.read(b); err != nil { // read required bytes amount counting taken bytes internally
		return 0, err
	}
	b1 := int(b[0])
	b2 := int(b[1])
	value := (0xff&b1)<<8 | (0xff & b2)
	return value, nil
}

func (d *decoder) readByte() (int, error) {
	b := allocateBytes(1)
	if err := d.read(b); err != nil { // read required bytes amount counting taken bytes internally
		return 0, err
	}

	byte1 := int(b[0])

	value := 0xff & byte1
	return value, nil
}

func (d *decoder) readBytes(size int) ([]byte, error) {
	b := allocateBytes(size)
	if err := d.read(b); err != nil { // read required bytes amount counting taken bytes internally
		return nil, err
	}

	return b, nil
}

func intSign(power int) int { /* "sign" power */
	num := -1 /* sign return value */

	if power == 0 {
		return 1
	}

	for cnt := 1; cnt < power; cnt++ {
		num *= -1
	}

	return num
}
