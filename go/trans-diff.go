// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"strings"
	"unsafe"

	"github.com/biogo/biogo/feat/genome/mouse/mm10"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/boom"
	"github.com/biogo/store/interval"
)

var (
	in,
	annot string
	classes set

	thresh float64

	pretty bool

	binLength int
	minLength int
	maxLength int

	minId  int
	minQ   int
	minAvQ float64
	mapQ   int
	mapQb  byte

	format string
)

const readLength = 50

const (
	all = iota
	primary
	secondary
)

type set []string

func (s *set) String() string {
	if len(*s) == 0 {
		return `""`
	}
	return strings.Join(*s, ",")
}
func (s *set) Set(value string) error {
	*s = append(*s, strings.Split(value, ",")...)
	if len(*s) == 0 {
		return errors.New("empty set")
	}
	return nil
}

func annotOK(annot string, classes []string) bool {
	if annot == "" && len(classes) == 0 {
		return true
	}
	return annot != "" && len(classes) != 0
}

func init() {
	flag.StringVar(&in, "in", "", "BAM file to be processed.")
	flag.StringVar(&annot, "annot", "", "file name of a GFF file containing annotations.")
	flag.Float64Var(&thresh, "thresh", 1, "log score threshold for inclusion of feature.")
	flag.Var(&classes, "class", "comma separated set of annotation classes to analyse.")
	flag.BoolVar(&pretty, "pretty", true, "outfile JSON data indented.")
	flag.IntVar(&minLength, "min", 20, "minimum length read considered.")
	flag.IntVar(&maxLength, "max", 35, "maximum length read considered.")
	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for mapped bases.")
	flag.IntVar(&minQ, "minQ", 20, "minimum per-base sequence quality.")
	flag.Float64Var(&minAvQ, "minAvQ", 30, "minimum average per-base sequence quality.")
	flag.IntVar(&mapQ, "mapQ", 0, "minimum mapping quality [0, 255).")
	flag.IntVar(&binLength, "bin", 1e7, "bin length.")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	mapQb = byte(mapQ)
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if in == "" || !annotOK(annot, classes) || mapQ < 0 || mapQ > 254 {
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

var index = map[string]int{}

func init() {
	for i, c := range mm10.Chromosomes {
		index[strings.ToLower(c.Chr)] = i
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

type intGff struct {
	*gff.Feature
}

func (f intGff) Range() interval.IntRange { return interval.IntRange{f.Start(), f.End()} }
func (f intGff) Overlap(b interval.IntRange) bool {
	return f.Feature.FeatEnd > b.Start && f.Feature.FeatStart < b.End
}
func (f intGff) ID() uintptr { return *(*uintptr)(unsafe.Pointer(f.Feature)) }

func filterFeats(annot string, classes, names []string, thresh float64) ([]interval.IntTree, error) {
	ntab := make(map[string]int, len(names))
	for i, n := range names {
		ntab[n] = i
	}

	ts := make([]interval.IntTree, len(names))

	cm := make(map[string]struct{})
	for _, c := range classes {
		cm[c] = struct{}{}
	}

	f, err := os.Open(annot)
	if err != nil {
		return nil, err
	}
	fs := featio.NewScanner(gff.NewReader(f))
	for fs.Next() {
		f := fs.Feat().(*gff.Feature)
		if f.FeatScore == nil || math.Exp(*f.FeatScore) < thresh {
			continue
		}
		var class string
		if att := f.FeatAttributes.Get(f.Feature); att != "" {
			// This gets the repeat attributes only.
			class = f.Feature + "/" + strings.Fields(att)[1]
		} else {
			class = f.Feature
		}
		if _, ok := cm[class]; !ok {
			last := strings.LastIndex(class, "/")
			if last == strings.Index(class, "/") {
				continue
			}
			if _, ok := cm[class[:last]]; !ok {
				continue
			}
		}
		if chr, ok := ntab[f.SeqName]; ok {
			ts[chr].Insert(intGff{f}, true)
		}
	}
	if err := fs.Error(); err != nil {
		return nil, err
	}
	for i := range ts {
		ts[i].AdjustRanges()
	}

	return ts, nil
}

type intBam struct {
	*boom.Record
}

func (f intBam) Range() interval.IntRange { return interval.IntRange{f.Start(), f.End()} }
func (f intBam) Overlap(b interval.IntRange) bool {
	return f.Start()+len(f.Seq()) > b.Start && f.Start() < b.End
}
func (f intBam) ID() uintptr { return *(*uintptr)(unsafe.Pointer(f.Record)) }

func main() {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer bf.Close()
	names := bf.RefNames()

	var classFilt []interval.IntTree
	if annot != "" {
		classFilt, err = filterFeats(annot, classes, names, thresh)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
	}

	var reads, totals int

	for {
		r, _, err := bf.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		if r.Flags()&boom.Unmapped == 0 {
			if qualOk(r, minId, minQ, minAvQ) {
				totals++
				if classFilt != nil && len(classFilt[r.RefID()].Get(intBam{r})) == 0 {
					continue
				}
				if s := r.Score(); s >= mapQb && s != 0xff {
					reads++
				}
			}
		}
	}

	fmt.Println(in, reads, totals, float64(reads)/float64(totals))
}
