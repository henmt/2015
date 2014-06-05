// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// length-heat generates a json description of piRNA expression levels
// by genomic location and piRNA length.
//
// A number of parameterised options are provided that allow tailoring of the analysis:
//
//  - piRNA type filtering;
//  - mapping quality filtering;
//  - piRNA deduplication by denesting reads;
//  - genomic bin size adjustment.
//
// Approach
//
// BAM alignments are read and filtered on sequence and mapping quality. Alignments that
// pass these filters are then filtered optionally on piRNA 1° or 2° status.
//
// Alignments are then counted, recording their length and position and unique 5' ends are
// counted as a proxy for piRNA family since piRNAs appear to truncate primarily from the 3'
// end. If denesting is requested, 5' ends of alignments are only considered if the alignment
// is not fully encompassed by another alignment.
//
// Denesting
//
// Denesting is performed by keeping a record of all unique alignment intervals in an interval
// tree and then after reading all alignments, unique 5' ends are identified as intervals with
// no containing interval.
package main

import (
	"code.google.com/p/biogo.boom"
	"code.google.com/p/biogo.store/interval"
	"code.google.com/p/biogo/feat"
	"code.google.com/p/biogo/feat/genome/mouse/mm10"

	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
)

var (
	in, out string
	pretty  bool

	filter int
	strict bool

	binLength int
	minLength int
	maxLength int

	minId  int
	minQ   int
	minAvQ float64
	mapQ   int
	mapQb  byte

	denest bool
)

const (
	all = iota
	primary
	secondary
)

func init() {
	flag.StringVar(&in, "in", "", "file name of a BAM file to be processed.")
	flag.StringVar(&out, "out", "", "outfile name.")
	flag.BoolVar(&pretty, "pretty", true, "outfile JSON data indented.")
	flag.BoolVar(&denest, "denest", false, "only consider denested reads for support count.")
	flag.IntVar(&minLength, "min", 20, "minimum length read considered.")
	flag.IntVar(&maxLength, "max", 35, "maximum length read considered.")
	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for non-clipped bases.")
	flag.IntVar(&minQ, "minQ", 20, "minimum per-base sequence quality.")
	flag.IntVar(&filter, "f", 0, "filter on piwi type 0: no filter, 1: primary, 2: secondary.")
	flag.BoolVar(&strict, "strict", false, "filter rejects ambiguous reads.")
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
				if t == boom.CigarMatch || t == boom.CigarEqual || t == boom.CigarSoftClipped {
					if t != boom.CigarSoftClipped {
						match++
					}
					mQ += int(q)
				}
			}
		}
	}
	match -= edit

	return match*100 >= minId*l && mQ >= int(minAvQ*float64(l))
}

type location struct {
	rid int
	pos int
}

type mappings struct {
	reads map[int]int
	kinds map[int]struct{}
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

type readSetElement struct {
	refid, start, length int
	reverse              boom.Flags
}

type read struct {
	*boom.Record
	id uintptr
}

func (r read) Range() interval.IntRange { return interval.IntRange{r.Start(), r.End()} }
func (r read) Overlap(b interval.IntRange) bool {
	// Half-open interval indexing.
	return r.End() > b.Start && r.Start() < b.End
}
func (r read) ID() uintptr { return r.id }
func (r read) End() int    { return r.Start() + len(r.Seq()) }

type containedRead struct {
	read
}

func (r containedRead) Overlap(b interval.IntRange) bool {
	// Return whether r is contained within b such that r is shorter than b.
	return (r.Start() > b.Start && r.End() <= b.End) || (r.Start() >= b.Start && r.End() < b.End)
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

func main() {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer bf.Close()

	bd := make(map[location]mappings)

	readSet := make(map[readSetElement]struct{})
	ts := make(map[int][]interval.IntTree)
loop:
	for id := uintptr(0); ; id++ {
		var r *boom.Record
		r, _, err = bf.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			break
		}
		if r.Flags()&boom.Unmapped == 0 {
			if qualOk(r, minId, minQ, minAvQ) {
				seq := r.Seq()
				if s, l := r.Score(), len(seq); s >= mapQb && s != 0xff && minLength <= l && l <= maxLength {
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
						panic("illegal filter value")
					}
					var sc mappings
					loc := location{rid: r.RefID(), pos: r.Start() / binLength}
					sc, ok := bd[loc]
					if !ok {
						sc.reads = make(map[int]int)
						sc.kinds = make(map[int]struct{})
					}
					sc.reads[l]++
					if denest {
						t, ok := ts[r.RefID()]
						if !ok {
							t = make([]interval.IntTree, 2)
							ts[r.RefID()] = t
						}
						re := readSetElement{r.RefID(), r.Start(), len(r.Seq()), r.Flags() & boom.Reverse}
						if _, ok := readSet[re]; !ok {
							readSet[re] = struct{}{}
							t[(r.Flags()&boom.Reverse)>>4].Insert(read{Record: r, id: id}, true)
						}
					} else {
						// Sign gives us the capacity to distinguish strands except at 0. This
						// is safe because the only way a reverse sense alignment can collide with
						// this is if the alignment length is zero and the start is at zero; this
						// cannot happen.
						if r.Flags()&boom.Reverse == 0 {
							sc.kinds[r.Start()] = struct{}{}
						} else {
							sc.kinds[-(r.Start() + l)] = struct{}{}
						}
					}
					bd[loc] = sc
				}
			}
		}
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	if denest {
		for rid, t := range ts {
			for i := range t {
				t[i].AdjustRanges()
				t[i].Do(func(iv interval.IntInterface) (done bool) {
					if c := len(t[i].Get(containedRead{iv.(read)})); c == 0 {
						r := iv.Range()
						loc := location{rid: rid, pos: r.Start / binLength}
						sc, ok := bd[loc]
						if !ok {
							panic("internal inconsistency")
						}
						if i == 0 {
							sc.kinds[r.Start] = struct{}{}
						} else {
							sc.kinds[-r.End] = struct{}{}
						}
						bd[loc] = sc
					}
					return
				})
				ts[rid][i] = interval.IntTree{}
			}
		}
	}

	var rna []*feature
	names := bf.RefNames()
	for k, v := range bd {
		i, ok := index[strings.ToLower(names[k.rid])]
		if !ok {
			continue
		}
		c := mm10.Chromosomes[i]
		f := &feature{
			start:    k.pos * binLength,
			end:      min((k.pos+1)*binLength, c.Len()),
			chr:      c,
			scores:   make([]float64, maxLength-minLength+1),
			supports: len(v.kinds),
		}
		for l, sv := range v.reads {
			f.scores[l-minLength] = float64(sv) / float64(f.Len())
		}
		rna = append(rna, f)
	}

	err = writeJSON(out, rna, filter, pretty)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
}

func decorate(out, format string, filter int) string {
	switch filter {
	case all:
		return fmt.Sprintf("%s.%s", out, format)
	case primary:
		return fmt.Sprintf("%s-U1.%s", out, format)
	case secondary:
		return fmt.Sprintf("%s-A10.%s", out, format)
	default:
		panic("illegal filter")
	}
}

func writeJSON(out string, rna []*feature, filter int, pretty bool) error {
	jsf, err := os.Create(decorate(out, "json", filter))
	if err != nil {
		return err
	}
	defer jsf.Close()

	type ranged struct {
		Sample string `json:"sample"`

		Bin    int `json:"bin"`
		Filter int `json:"filter"`

		Min int `json:"min"`
		Max int `json:"max"`

		MinQ   int     `json:"min-qual"`
		MinAvQ float64 `json:"min-av-qual"`
		MinID  int     `json:"min-id"`
		MapQ   int     `json:"map-qual"`

		Features []*feature `json:"features"`
	}

	r := ranged{
		Sample:   path(in),
		Bin:      binLength,
		Filter:   filter,
		Min:      minLength,
		Max:      maxLength,
		MinQ:     minQ,
		MinAvQ:   minAvQ,
		MinID:    minId,
		MapQ:     mapQ,
		Features: rna,
	}

	if pretty {
		j, err := json.MarshalIndent(r, "", "  ")
		if err != nil {
			return err
		}
		_, err = jsf.Write(j)
		if err != nil {
			return err
		}
	} else {
		enc := json.NewEncoder(jsf)
		err = enc.Encode(r)
		if err != nil {
			return err
		}
	}

	return nil
}

func path(p string) string {
	p, _ = filepath.Abs(p)
	return p
}

type feature struct {
	start, end int
	name       string
	chr        feat.Feature
	scores     []float64
	supports   int
}

func (f *feature) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Chr      string    `json:"chr"`
		Start    int       `json:"start"`
		End      int       `json:"end"`
		Scores   []float64 `json:"scores"`
		Supports int       `json:"support"`
	}{
		Chr:      f.chr.Name(),
		Start:    f.start,
		End:      f.end,
		Scores:   f.scores,
		Supports: f.supports,
	})
}

func (f *feature) Start() int             { return f.start }
func (f *feature) End() int               { return f.end }
func (f *feature) Len() int               { return f.end - f.start }
func (f *feature) Name() string           { return f.name }
func (f *feature) Description() string    { return "alignment bin" }
func (f *feature) Location() feat.Feature { return f.chr }
func (f *feature) Scores() []float64      { return f.scores }
