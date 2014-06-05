// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// length-heat-annot-diff generates a json description of piRNA expression level differences
// by genomic location and piRNA length.
//
// A number of parameterised options are provided that allow tailoring of the analysis:
//
//  - piRNA type filtering;
//  - feature overlap filtering;
//  - mapping quality filtering;
//  - piRNA deduplication by denesting reads;
//  - genomic bin size adjustment.
//
// Approach
//
// Two BAM alignments are read and filtered on sequence and mapping quality. Alignments that
// pass these filters are then filtered optionally on piRNA 1° or 2° status.
//
// Alignments that pass these filters are then assessed for overlap with an optionally provided
// set of features from a GFF file and a set of feature classes to compare against. If these are
// given only alignments that overlap the specified classes of features are included in the
// analysis.
//
// Alignments are then counted, recording their length and position and unique 5' ends are
// counted as a proxy for piRNA family since piRNAs appear to truncate primarily from the 3'
// end. If denesting is requested, 5' ends of alignments are only considered if the alignment
// is not fully encompassed by another alignment.
//
// The difference between the two inputs' bin and length tallies are then calculated and stored.
// Unique 5' end counts are kept as individual data from each sample.
//
// Assumptions
//
// We are interested in the source of the piRNA rather than the targets of the piRNA in this
// instance.
//
// The source will match the read better than the target and that only a small number of
// mismatching reads will be reported (and further that if they mismatch too much they will
// be culled by % id filtering).
//
// In the case of repeats, some proportion of the piRNA source loci will be surrounded by enough
// target-like sequence that it will be recognised by RepeatMasker and so will be annotated.
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
	"code.google.com/p/biogo/feat/genome"
	"code.google.com/p/biogo/feat/genome/mouse/mm10"
	"code.google.com/p/biogo/io/featio"
	"code.google.com/p/biogo/io/featio/gff"

	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"unsafe"
)

var (
	in pair
	annot,
	ref,
	out string
	classes set

	pretty bool

	filter int

	binLength int
	minLength int
	maxLength int

	minId  int
	minQ   int
	minAvQ float64
	mapQ   int
	mapQb  byte

	denest bool

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

type pair [2]string

func (p *pair) String() string {
	return fmt.Sprintf("%q,%q", p[0], p[1])
}
func (p *pair) Set(value string) error {
	c := strings.Split(value, ",")
	switch len(c) {
	case 0:
		return errors.New("empty pair")
	case 2:
		copy((*p)[:], c)
		return nil
	default:
		return fmt.Errorf("unexpected number of elements: got %d expected 2", len(c))
	}
}

func annotOK(annot string, classes []string) bool {
	if annot == "" && len(classes) == 0 {
		return true
	}
	return annot != "" && len(classes) != 0
}

func init() {
	flag.Var(&in, "in", "comma separated pair of BAM files to be processed.")
	flag.StringVar(&ref, "ref", "", "fasta file of the genome to be processed.")
	flag.StringVar(&annot, "annot", "", "file name of a GFF file containing annotations.")
	flag.Var(&classes, "class", "comma separated set of annotation classes to analyse.")
	flag.StringVar(&out, "out", "", "outfile name.")
	flag.BoolVar(&pretty, "pretty", true, "outfile JSON data indented.")
	flag.BoolVar(&denest, "denest", false, "only consider denested reads for support count.")
	flag.IntVar(&minLength, "min", 20, "minimum length read considered.")
	flag.IntVar(&maxLength, "max", 35, "maximum length read considered.")
	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for mapped bases.")
	flag.IntVar(&minQ, "minQ", 20, "minimum per-base sequence quality.")
	flag.IntVar(&filter, "f", 0, "filter on piwi type 0: no filter, 1: primary, 2: secondary.")
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
	if in[0] == "" || in[1] == "" || out == "" || !annotOK(annot, classes) || mapQ < 0 || mapQ > 254 {
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

func checkNames(inPair pair) ([]string, error) {
	var names []string
	for _, in := range inPair {
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

func filterFeats(annot string, classes, names []string) ([]interval.IntTree, error) {
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

type location struct {
	rid int
	bin int
}

type bins struct {
	totals   [2]int
	mappings map[location]mappings
}

type mappings struct {
	reads [][2]int
	kinds [2]map[int]struct{}
}

func (b *bins) merge(a *bins) *bins {
	if b == a {
		return b
	}
	if b == nil {
		return a
	}
	for l, al := range a.mappings {
		if _, ok := b.mappings[l]; !ok {
			b.mappings[l] = al
		} else {
			bl := b.mappings[l]

			for i, v := range al.reads {
				bl.reads[i][0] += v[0]
				bl.reads[i][1] += v[1]
			}

			for i, v := range al.kinds {
				if bl.kinds[i] == nil {
					bl.kinds[i] = v
				} else {
					for k := range v {
						bl.kinds[i][k] = struct{}{}
					}
				}
			}

			b.mappings[l] = bl
		}
	}

	b.totals[0] += a.totals[0]
	b.totals[1] += a.totals[1]

	return b
}

func rnaFeats(inPair pair, names []string, classFilt []interval.IntTree, minLength, maxLength, minId, minQ int, mapQb byte, minAvQ float64) (sf []*feature, totals [2]int, err error) {
	defer func() {
		if r := recover(); r != nil {
			sf, err = nil, r.(error)
		}
	}()

	var (
		bd *bins
		rc = make(chan *bins)
		wg sync.WaitGroup
	)
	go func() {
		for i, in := range inPair {
			wg.Add(1)
			go func(id int, in string) {
				defer wg.Done()
				var b *bins
				b, err = readBam(id, in, classFilt, minLength, maxLength, minId, minQ, mapQb, minAvQ)
				if err != nil {
					panic(err)
				}
				rc <- b
			}(i, in)
		}
		wg.Wait()
		close(rc)
	}()

	for b := range rc {
		bd = bd.merge(b)
	}

	for loc, scores := range bd.mappings {
		i, ok := index[strings.ToLower(names[loc.rid])]
		if !ok {
			continue
		}
		c := mm10.Chromosomes[i]
		f := &feature{
			start: loc.bin * binLength,
			end:   min((loc.bin+1)*binLength, c.Len()),
			typ:   "delta",
			chr:   c,
			counts: [2][]int{
				make([]int, maxLength-minLength+1),
				make([]int, maxLength-minLength+1),
			},
			supports: [2]int{len(scores.kinds[0]), len(scores.kinds[1])},
		}
		for i, pv := range scores.reads {
			f.counts[0][i] = pv[0]
			f.counts[1][i] = pv[1]
		}
		sf = append(sf, f)
	}

	// Fill any holes due to absence of reads.
	for rid, name := range names {
		i, ok := index[strings.ToLower(name)]
		if !ok {
			continue
		}
		c := mm10.Chromosomes[i]
		for bin := c.Start(); bin*binLength < c.End(); bin++ {
			if _, ok := bd.mappings[location{rid: rid, bin: bin}]; !ok {
				sf = append(sf, &feature{
					start: bin * binLength,
					end:   min((bin+1)*binLength, c.Len()),
					typ:   "missing",
					chr:   c,
					counts: [2][]int{
						make([]int, maxLength-minLength+1),
						make([]int, maxLength-minLength+1),
					},
				})
			}
		}
	}

	return sf, bd.totals, nil
}

func (f intBam) Range() interval.IntRange { return interval.IntRange{f.Start(), f.End()} }
func (f intBam) Overlap(b interval.IntRange) bool {
	return f.Start()+len(f.Seq()) > b.Start && f.Start() < b.End
}
func (f intBam) ID() uintptr { return *(*uintptr)(unsafe.Pointer(f.Record)) }

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

func readBam(id int, in string, classFilt []interval.IntTree, minLength, maxLength, minId, minQ int, mapQb byte, minAvQ float64) (*bins, error) {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		return nil, err
	}
	defer bf.Close()

	bd := &bins{mappings: make(map[location]mappings)}

	readSet := [2]map[readSetElement]struct{}{
		make(map[readSetElement]struct{}),
		make(map[readSetElement]struct{}),
	}
	ts := make(map[int][][2]interval.IntTree)
loop:
	for uid := uintptr(0); ; uid++ {
		r, _, err := bf.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if r.Flags()&boom.Unmapped == 0 {
			if qualOk(r, minId, minQ, minAvQ) {
				bd.totals[id]++
				seq := r.Seq()
				if s, l := r.Score(), len(seq); s >= mapQb && s != 0xff && minLength <= l && l <= maxLength {
					switch filter {
					case all:
					case primary:
						if l < 1 {
							continue loop
						}
						if r.Flags()&boom.Reverse == 0 {
							if seq[0]|' ' != 't' {
								continue loop
							}
						} else {
							if seq[l-1]|' ' != 'a' {
								continue loop
							}
						}
					case secondary:
						if l < 10 {
							continue loop
						}
						if r.Flags()&boom.Reverse == 0 {
							if seq[9]|' ' != 'a' {
								continue loop
							}
						} else {
							if seq[l-10]|' ' != 't' {
								continue loop
							}
						}
					default:
						panic("illegal filter value")
					}

					if classFilt != nil && len(classFilt[r.RefID()].Get(intBam{r})) == 0 {
						continue loop
					}

					loc := location{rid: r.RefID(), bin: r.Start() / binLength}
					sc, ok := bd.mappings[loc]
					if !ok {
						sc.reads = make([][2]int, maxLength-minLength+1)
						for i := range sc.kinds {
							sc.kinds[i] = make(map[int]struct{})
						}
					}
					sc.reads[l-minLength][id]++
					t, ok := ts[r.RefID()]
					if !ok {
						t = make([][2]interval.IntTree, 2)
						ts[r.RefID()] = t
					}
					if denest {
						re := readSetElement{r.RefID(), r.Start(), len(r.Seq()), r.Flags() & boom.Reverse}
						if _, ok := readSet[id][re]; !ok {
							readSet[id][re] = struct{}{}
							t[id][(r.Flags()&boom.Reverse)>>4].Insert(read{Record: r, id: uid}, true)
						}
					} else {
						// Sign gives us the capacity to distinguish strands except at 0. This
						// is safe because the only way a reverse sense alignment can collide with
						// this is if the alignment length is zero and the start is at zero; this
						// cannot happen.
						if r.Flags()&boom.Reverse == 0 {
							sc.kinds[id][r.Start()] = struct{}{}
						} else {
							sc.kinds[id][-(r.Start() + l)] = struct{}{}
						}
					}
					bd.mappings[loc] = sc
				}
			}
		}
	}

	if denest {
		for rid, t := range ts {
			for id := range t {
				for strand := range t[id] {
					t[id][strand].AdjustRanges()
					t[id][strand].Do(func(iv interval.IntInterface) (done bool) {
						if c := len(t[id][strand].Get(containedRead{iv.(read)})); c == 0 {
							r := iv.Range()
							loc := location{rid: rid, bin: r.Start / binLength}
							sc, ok := bd.mappings[loc]
							if !ok {
								panic("internal inconsistency")
							}
							if strand == 0 {
								sc.kinds[id][r.Start] = struct{}{}
							} else {
								sc.kinds[id][-r.End] = struct{}{}
							}
							bd.mappings[loc] = sc
						}
						return
					})
				}
			}
		}
	}

	return bd, nil
}

func main() {
	names, err := checkNames(in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	var classFilt []interval.IntTree
	if annot != "" {
		classFilt, err = filterFeats(annot, classes, names)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
	}

	rna, totals, err := rnaFeats(in, names, classFilt, minLength, maxLength, minId, minQ, mapQb, minAvQ)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	err = writeJSON(out, rna, totals, filter, pretty)
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

func writeJSON(out string, rna []*feature, totals [2]int, filter int, pretty bool) error {
	jsf, err := os.Create(decorate(out, "json", filter))
	if err != nil {
		return err
	}
	defer jsf.Close()

	type ranged struct {
		Pair [2]string `json:"pair"`

		Bin     int      `json:"bin"`
		Classes []string `json:"classes"`
		Filter  int      `json:"filter"`

		Min int `json:"min"`
		Max int `json:"max"`

		MinQ   int     `json:"min-qual"`
		MinAvQ float64 `json:"min-av-qual"`
		MinID  int     `json:"min-id"`
		MapQ   int     `json:"map-qual"`

		Totals   [2]int     `json:"totals"`
		Features []*feature `json:"features"`
	}

	r := ranged{
		Pair:     paths(in),
		Bin:      binLength,
		Classes:  classes,
		Filter:   filter,
		Min:      minLength,
		Max:      maxLength,
		MinQ:     minQ,
		MinAvQ:   minAvQ,
		MinID:    minId,
		MapQ:     mapQ,
		Totals:   totals,
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

func paths(p pair) [2]string {
	var a [2]string
	for i := range p {
		a[i], _ = filepath.Abs(p[i])
	}
	return a
}

type feature struct {
	chr *genome.Chromosome
	start,
	end int

	typ string

	counts   [2][]int
	supports [2]int
}

func (f *feature) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Chr      string   `json:"chr"`
		Start    int      `json:"start"`
		End      int      `json:"end"`
		Type     string   `json:"type"`
		Counts   [2][]int `json:"counts"`
		Supports [2]int   `json:"support"`
	}{
		Chr:      f.chr.Name(),
		Start:    f.start,
		End:      f.end,
		Type:     f.typ,
		Counts:   f.counts,
		Supports: f.supports,
	})
}
