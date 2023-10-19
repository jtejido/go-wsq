package main

import (
	"fmt"
	"image/jpeg"
	"log"
	"os"
	"wsq"
)

// simple decode->encode->decode then convert to jpeg test
func main() {
	f, err := os.Open("sample.wsq")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	img, err := wsq.Decode(f)
	if err != nil {
		log.Fatal(err)
	}

	fo, err := os.Create("encoded.wsq")
	if err != nil {
		log.Fatal(err)
	}

	if err := wsq.Encode(fo, img, &wsq.Options{
		Comments: []string{"test"},
	}); err != nil {
		fmt.Println("err")
		log.Fatal(err)
	}
	fo.Close()
	f2, err := os.Open("encoded.wsq")
	if err != nil {
		log.Fatal(err)
	}
	defer f2.Close()

	img2, err := wsq.Decode(f2)

	if err != nil {
		log.Fatal(err)
	}

	ff, err := os.Create("img.jpeg")
	if err != nil {
		panic(err)
	}
	defer ff.Close()
	if err := jpeg.Encode(ff, img2, nil); err != nil {
		log.Fatal(err)
	}
}
