// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/feat/genome/mouse/mm10"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/boom"
	"github.com/biogo/store/interval"
)

const (
	all = iota
	primary
	secondary
)

var (
	annot string
	reads,
	classes set
	filter int
	strict bool
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

func init() {
	flag.Var(&reads, "reads", "comma separated set of BAM file to be processed.")
	flag.StringVar(&annot, "annot", "", "file name of a GFF file containing annotations.")
	flag.Var(&classes, "class", "comma separated set of annotation classes to analyse.")
	flag.IntVar(&filter, "f", 0, "filter on piwi type 0: no filter, 1: primary, 2: secondary.")
	flag.BoolVar(&strict, "strict", false, "filter rejects ambiguous reads.")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if len(reads) == 0 || annot == "" {
		flag.Usage()
		os.Exit(1)
	}
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
	uintptr
}

func (f intGff) Range() interval.IntRange { return interval.IntRange{f.Start(), f.End()} }
func (f intGff) Overlap(b interval.IntRange) bool {
	return f.Feature.FeatEnd > b.Start && f.Feature.FeatStart < b.End
}
func (f intGff) ID() uintptr { return f.uintptr }

func annotFeats(annot string, classes, names []string) ([]interval.IntTree, error) {
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
	for id := uintptr(0); fs.Next(); id++ {
		f := fs.Feat().(*gff.Feature)
		var class string
		att := f.FeatAttributes.Get("repeat")
		if att == "" {
			// Ignore non-repeat features.
			continue
		}
		repeatFields := strings.Fields(att)
		class = f.Feature + "/" + repeatFields[1]
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
			ts[chr].Insert(intGff{f, id}, true)
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
func (f intBam) ID() uintptr { return 0 }

type gffFeatures []*gff.Feature

func (f gffFeatures) Len() int           { return len(f) }
func (f gffFeatures) Less(i, j int) bool { return *f[i].FeatScore > *f[j].FeatScore }
func (f gffFeatures) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

func checkNames(files set) ([]string, error) {
	var names []string
	for _, in := range files {
		bf, err := boom.OpenBAM(in)
		if err != nil {
			return nil, err
		}
		if names != nil {
			for i, n := range bf.RefNames() {
				if names[i] != n {
					return nil, errors.New("header mismatch")
				}
			}
		}
		names = bf.RefNames()
		bf.Close()
	}
	return names, nil
}

func isPrimary(r *boom.Record) bool {
	seq := r.Seq()
	if len(seq) < 1 {
		return false
	}
	if r.Flags()&boom.Reverse == 0 {
		return seq[0]|' ' == 't'
	}
	return seq[len(seq)-1]|' ' == 'a'
}

func isSecondary(r *boom.Record) bool {
	seq := r.Seq()
	if len(seq) < 10 {
		return false
	}
	if r.Flags()&boom.Reverse == 0 {
		return seq[9]|' ' == 'a'
	}
	return seq[len(seq)-10]|' ' == 't'
}

type vector [][2]int

func (v *vector) inc(i int, s seq.Strand) {
	if s == 0 {
		return
	}
	if s < 0 {
		s = 0
	}
	switch {
	case i < 0:
	case i < len(*v):
		(*v)[i][s]++
	case i == len(*v):
		var p [2]int
		p[s] = 1
		*v = append(*v, p)
	default:
		t := make([][2]int, i+1)
		copy(t, *v)
		t[i][s] = 1
		*v = t
	}
}

func main() {
	names, err := checkNames(reads)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	feats, err := annotFeats(annot, classes, names)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	vm := make(map[string]*vector)
	fm := make(map[string]string)

	for _, in := range reads {
		fmt.Fprintf(os.Stderr, "Reading %q\n", in)
		bf, err := boom.OpenBAM(in)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}

	loop:
		for {
			r, _, err := bf.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}

			switch filter {
			case all:
			case primary:
				if !isPrimary(r) {
					continue loop
				}
				if strict && isSecondary(r) {
					continue loop
				}
			case secondary:
				if !isSecondary(r) {
					continue loop
				}
				if strict && isPrimary(r) {
					continue loop
				}
			default:
				panic(fmt.Errorf("illegal filter %d", filter))
			}

			if r.Flags()&boom.Unmapped == 0 {
				feats[r.RefID()].DoMatching(func(iv interval.IntInterface) (done bool) {
					f := iv.(intGff)
					repeatFields := strings.Fields(f.FeatAttributes.Get("repeat"))
					if len(repeatFields) == 0 {
						return
					}
					from, err := strconv.Atoi(repeatFields[2])
					if err != nil {
						panic(err)
					}
					v, ok := vm[repeatFields[0]]
					if !ok {
						v = &vector{}
						vm[repeatFields[0]] = v
						fm[repeatFields[0]] = repeatFields[1]
					}
					for pos := r.Start(); pos < r.End(); pos++ {
						var loc int
						switch f.FeatStrand {
						case seq.Plus:
							loc = pos - f.FeatStart + from
						case seq.Minus:
							loc = f.FeatEnd - pos + from
						default:
							return
						}
						var st seq.Strand
						if r.Flags()&boom.Reverse != 0 {
							st = -f.FeatStrand
						} else {
							st = f.FeatStrand
						}
						v.inc(loc, st)
					}
					return
				}, intBam{r})
			}
		}

		bf.Close()
	}

	for typ, vec := range vm {
		for pos, val := range *vec {
			if val[0] != 0 {
				fmt.Printf("%s\t%s\tminus\t%d\t%d\n", fm[typ], typ, pos, val[0])
			}
			if val[1] != 0 {
				fmt.Printf("%s\t%s\tplus\t%d\t%d\n", fm[typ], typ, pos, val[1])
			}
		}
	}
}
