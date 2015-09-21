// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/biogo/boom"
)

var (
	in, out string

	minId  int
	minQ   int
	minAvQ float64
	mapQ   int
	mapQb  byte
)

func init() {
	flag.StringVar(&in, "in", "", "file name of a BAM file to filter.")
	flag.StringVar(&out, "out", "", "outfile name.")
	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for mapped bases.")
	flag.IntVar(&minQ, "minQ", 20, "minimum per-base sequence quality.")
	flag.Float64Var(&minAvQ, "minAvQ", 30, "minimum average per-base sequence quality.")
	flag.IntVar(&mapQ, "mapQ", 0, "minimum mapping quality [0, 255).")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	mapQb = byte(mapQ)
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if in == "" || out == "" || mapQ < 0 || mapQ > 254 {
		flag.Usage()
		os.Exit(1)
	}
}

func qualOk(r *boom.Record, minId, minQ int, minAvQ float64) (ok bool) {
	var (
		off, l int
		match  int
		mQ     int
		edit   int
	)

	for _, t := range r.Tags() {
		if t.Tag() == [2]byte{'N', 'M'} {
			switch e := t.Value().(type) {
			case byte:
				edit = int(e)
			case uint16:
				edit = int(e)
			case uint32:
				edit = int(e)
			default:
				edit = 0
			}
		}
	}
	cigar := r.Cigar()
	qual := r.Quality()
	for _, c := range cigar {
		t := c.Type()
		if t == boom.CigarMatch || t == boom.CigarInsertion || t == boom.CigarSoftClipped || t == boom.CigarEqual || t == boom.CigarMismatch {
			off = l
			l += c.Len()
			for _, q := range qual[off:l] {
				if int(q) < minQ {
					return false
				}
				if t == boom.CigarMatch || t == boom.CigarEqual {
					match++
					mQ += int(q)
				}
			}
		}
	}
	match -= edit

	ok = match*100 >= minId*l && mQ >= int(minAvQ*float64(l))

	return
}

func filterBam(in, out string, minId, minQ int, mapQb byte, minAvQ float64) error {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		return err
	}
	defer bf.Close()

	bo, err := boom.CreateBAM(out, bf.Header(), true)
	if err != nil {
		return err
	}
	defer bo.Close()

	for {
		r, _, err := bf.Read()
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}
		if r.Score() >= mapQb && qualOk(r, minId, minQ, minAvQ) {
			_, err = bo.Write(r)
			if err != nil {
				return err
			}
		}
	}
}

func main() {
	err := filterBam(in, out, minId, minQ, mapQb, minAvQ)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
