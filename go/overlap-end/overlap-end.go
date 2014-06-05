// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// overlap-end generates a csv of values and associated plot describing truncation of piRNA
// alignments from two non-overlapping size pools
//
// A number of parameterised options are provided that allow tailoring of the analysis:
//
//  - pool alignment lengths (short and long);
//  - piRNA type filtering;
//  - mapping quality filtering;
//  - long pool piRNA deduplication by denesting reads; and
//  - short alignment containment.
//
// Approach
//
// Read alignments from one or two BAM files - if two files are specified one is nominally
// wild-type and the other nominally mutant, long alignments (default 28-32bp) are taken
// from the wild-type data and short alignments (default 23-27bp) are taken from the mutant
// data. These two data sets are then compared as follows.
//
// One instance of each unique alignment is kept from the longer pool, optionally denesting the
// unique long pool alignments. The initial uniqueness test is done by start/length identity
// (map assignment).
//
// For each alignment from the shorter pool find all long alignments that overlap their map
// location by at least one base in the non-contained case or all long alignments that completely
// overlap it in the contained case.
//
// Record the relative end positions between the short alignment and the found long alignments
// (transformed so that lengthening would be a positive distance and shortening would be a negative
// distance).
//
// Denesting
//
// Denesting is performed by keeping a record of all unique alignment intervals in an interval
// tree and then after reading all alignments, canonical long piRNA alignments are identified
// as intervals with no containing interval.
//
// Containment
//
// Containment of short reads is intended to reduce the linear scoring behaviour in response to
// long pool read truncations; a long pool with a number of truncated forms that overlap with
// short alignments would otherwise multiply out the counts of truncations. It is also intended
// to abolish a 'negative' truncation effect where short alignments are determined to be
// negatively truncated from a long alignment overlaps, but does not completely contain the short
// alignment.
//
// When contain is requested as an option (the default), long alignments are only considered for
// comparison to a short alignment when they completely overlap the short alignment and are longer
// than the short alignment.
package main

import (
	"code.google.com/p/biogo.boom"
	"code.google.com/p/biogo.store/interval"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/plotutil"
	"code.google.com/p/plotinum/vg"
	"code.google.com/p/plotinum/vg/vgsvg"

	"errors"
	"flag"
	"fmt"
	"image/color"
	"io"
	"os"
	"path/filepath"
	"strings"
)

const (
	all = iota
	primary
	secondary
)

var (
	filter  int
	strict  bool
	care    bool
	denest  bool
	contain bool

	pairs pair
	out   string

	longMinLength int
	longMaxLength int

	shortMinLength int
	shortMaxLength int

	maxLength int

	minId  int
	minQ   int
	minAvQ float64
	mapQ   int
	mapQb  byte
)

type pair [][2]string

func (p *pair) String() string {
	return fmt.Sprint(*p)
}
func (p *pair) Set(value string) error {
	var (
		pv [2]string
		ps = pv[:0]
	)
	for i, e := range strings.Split(value, ",") {
		if i > 1 {
			return errors.New("pair: too many parts")
		}
		if e == "" {
			return errors.New("pair: empty value")
		}
		ps = append(ps, e)
	}
	switch len(ps) {
	case 0:
		return errors.New("pair: no value")
	case 1:
		ps = append(ps, ps[0])
	}
	*p = append(*p, pv)
	return nil
}

const readLength = 50

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func init() {
	flag.Var(&pairs, "pair", "either a comma-separated pair of BAM files (long first) to be\n\tprocessed, or a single BAM file to be used for long and short.\n\t(may be invoked multiple times.)")
	flag.StringVar(&out, "out", "", "base name for output files.")

	flag.IntVar(&shortMinLength, "shortmin", 23, "minimum length short read considered.")
	flag.IntVar(&shortMaxLength, "shortmax", 27, "maximum length short read considered.")

	flag.IntVar(&longMinLength, "longmin", 28, "minimum length long read considered.")
	flag.IntVar(&longMaxLength, "longmax", 32, "maximum length long read considered.")

	flag.IntVar(&filter, "f", 0, "filter on piwi type 0: no filter, 1: primary, 2: secondary.")
	flag.BoolVar(&strict, "strict", false, "filter rejects ambiguous reads.")
	flag.BoolVar(&care, "care", true, "care whether the short reads also satisfy filter.")
	flag.BoolVar(&denest, "denest", false, "remove long reads that are nested within another long read.")
	flag.BoolVar(&contain, "contain", false, "only consider long reads completely containing a short query.")

	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for read sequence bases.")
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
	if len(pairs) == 0 || out == "" || mapQ < 0 || mapQ > 254 {
		flag.Usage()
		os.Exit(1)
	}
	maxLength = max(shortMaxLength, longMaxLength)
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

type readSetElement struct {
	refid, start, length int
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

func longTrees(long string) ([][2]interval.IntTree, error) {
	bf, err := boom.OpenBAM(long)
	if err != nil {
		return nil, fmt.Errorf("%v: %v", err, long)
	}
	defer bf.Close()

	ts := make([][2]interval.IntTree, bf.Targets())
	readSet := make(map[readSetElement]struct{})
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
			re := readSetElement{r.RefID(), r.Start(), len(r.Seq())}
			if !(longMinLength <= re.length && re.length <= longMaxLength) {
				continue
			}
			if _, ok := readSet[re]; ok {
				continue
			}
			readSet[re] = struct{}{}

			if qualOk(r, minId, minQ, minAvQ) {
				if s := r.Score(); s >= mapQb && s != 0xff {
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
						return nil, fmt.Errorf("illegal filter %d", filter)
					}

					ts[r.RefID()][(r.Flags()&boom.Reverse)>>4].Insert(read{Record: r, id: id}, true)
				}
			}
		}
	}
	for _, t := range ts {
		for i := range t {
			t[i].AdjustRanges()

			if denest {
				var rm []interval.IntInterface
				t[i].Do(func(iv interval.IntInterface) (done bool) {
					if len(t[i].Get(containedRead{iv.(read)})) > 0 {
						rm = append(rm, iv)
					}
					return
				})

				for _, iv := range rm {
					t[i].Delete(iv, true)
				}
				t[i].AdjustRanges()
			}
		}
	}

	return ts, err
}

type set struct {
	long     string
	short    string
	fiveEnd  []int
	threeEnd []int
}

func searchForest(ts [][2]interval.IntTree, contain bool, short string) (set, error) {
	bf, err := boom.OpenBAM(short)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return set{}, fmt.Errorf("%v: %v", err, short)
	}
	defer bf.Close()

	results := set{
		fiveEnd:  make([]int, 2*maxLength),
		threeEnd: make([]int, 2*maxLength),
	}
loop:
	for {
		var r *boom.Record
		r, _, err = bf.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			break
		}
		if r.Flags()&boom.Unmapped == 0 {
			if seq := r.Seq(); !(shortMinLength <= len(seq) && len(seq) <= shortMaxLength) {
				continue
			}
			if qualOk(r, minId, minQ, minAvQ) {
				if s := r.Score(); s >= mapQb && s != 0xff {
					if care {
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
							return set{}, fmt.Errorf("illegal filter %d", filter)
						}
					}

					var longs []interval.IntInterface
					if contain {
						// Using containedRead is a little lazy, but we know the short read is shorter
						// than all the reads in the long tree, so we can use this type.
						longs = ts[r.RefID()][(r.Flags()&boom.Reverse)>>4].Get(containedRead{read{Record: r}})
					} else {
						longs = ts[r.RefID()][(r.Flags()&boom.Reverse)>>4].Get(read{Record: r})
					}
					for _, long := range longs {
						longr := long.(read)

						// These distances result in a shift to the left for shorter short reads.
						if r.Flags()&boom.Reverse == 0 {
							results.fiveEnd[longr.Start()-r.Start()+maxLength]++
							results.threeEnd[read{Record: r}.End()-longr.End()+maxLength]++
						} else {
							results.fiveEnd[read{Record: r}.End()-longr.End()+maxLength]++
							results.threeEnd[longr.Start()-r.Start()+maxLength]++
						}
					}
				}
			}
		}
	}

	return results, err
}

func main() {
	var data []set
	for _, p := range pairs {
		long, short := p[0], p[1]

		trees, err := longTrees(long)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}

		d, err := searchForest(trees, contain, short)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		d.long = filepath.Base(long)
		d.short = filepath.Base(short)

		data = append(data, d)

		if len(pairs) == 1 {
			if err := barchart(out, d); err != nil {
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}
		}
	}

	if err := csv(out, data); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	if len(data) > 1 {
		if err := boxplot(out, data); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
	}
}

func decorate(path, decoration string, filter int) string {
	switch filter {
	case all:
		return fmt.Sprintf("%s.%s", path, decoration)
	case primary:
		return fmt.Sprintf("%s-U1.%s", path, decoration)
	case secondary:
		return fmt.Sprintf("%s-A10.%s", path, decoration)
	default:
		panic("illegal filter")
	}
}

func csv(path string, data []set) error {
	if len(data) == 0 {
		return errors.New("no data")
	}

	f, err := os.Create(decorate(path, "csv", filter))
	if err != nil {
		return err
	}
	defer f.Close()

	_, err = fmt.Fprint(f, "End,Long,Short")
	if err != nil {
		return err
	}
	for i := range data[0].fiveEnd {
		_, err = fmt.Fprintf(f, ",%d", i-maxLength)
		if err != nil {
			return err
		}
	}
	_, err = fmt.Fprintln(f)
	if err != nil {
		return err
	}

	for _, e := range data {
		_, err = fmt.Fprintf(f, "FivePrime,%s,%s", e.long, e.short)
		if err != nil {
			return err
		}
		for _, v := range e.fiveEnd {
			_, err = fmt.Fprintf(f, ",%d", v)
			if err != nil {
				return err
			}
		}
		_, err = fmt.Fprintf(f, "\nThreePrime,%s,%s", e.long, e.short)
		if err != nil {
			return err
		}
		for _, v := range e.threeEnd {
			_, err = fmt.Fprintf(f, ",%d", v)
			if err != nil {
				return err
			}
		}
		_, err = fmt.Fprintln(f)
	}

	return err
}

type normalised struct {
	vals []int
	sum  float64
}

func (n *normalised) Len() int {
	n.sum = 0
	for _, v := range n.vals {
		n.sum += float64(v)
	}
	return len(n.vals)
}
func (n *normalised) Value(i int) float64 { return float64(n.vals[i]) / n.sum }

var titles = [...]string{
	"Read end offsets - all small RNA",
	"Read end offsets - primary piRNA",
	"Read end offsets - secondary piRNA",
}

func barchart(path string, data set) error {
	font, err := vg.MakeFont("Helvetica", 10)
	if err != nil {
		return err
	}
	titleFont, err := vg.MakeFont("Helvetica", 12)
	if err != nil {
		return err
	}
	style := plot.TextStyle{Color: color.Gray{0}, Font: font}
	p, err := plot.New()
	if err != nil {
		return err
	}
	p.Title.Text = titles[filter]
	p.Title.TextStyle = plot.TextStyle{Color: color.Gray{0}, Font: titleFont}
	p.X.Label.Text = "Length Offset"
	p.Y.Label.Text = "Relative Frequency"
	p.X.Label.TextStyle = style
	p.Y.Label.TextStyle = style
	p.X.Tick.Label = style
	p.Y.Tick.Label = style
	p.Legend.TextStyle = style

	barsFivePrime, err := plotter.NewBarChart(&normalised{vals: data.fiveEnd}, 1) // A non-zero width is required to prevent the creation failing.
	if err != nil {
		return err
	}
	barsFivePrime.LineStyle.Width = vg.Length(0)
	barsFivePrime.Color = plotutil.Color(0)

	barsThreePrime, err := plotter.NewBarChart(&normalised{vals: data.threeEnd}, 1) // A non-zero width is required to prevent the creation failing.
	if err != nil {
		return err
	}
	barsThreePrime.LineStyle.Width = vg.Length(0)
	barsThreePrime.Color = plotutil.Color(1)

	p.Add(barsFivePrime, barsThreePrime)
	p.Legend.Add("5'-end", barsFivePrime)
	p.Legend.Add("3'-end", barsThreePrime)
	p.Legend.Top = true
	p.NominalX(func() []string {
		n := make([]string, len(data.fiveEnd))
		for i := range data.fiveEnd {
			if v := i - maxLength; v%5 == 0 {
				n[i] = fmt.Sprint(v)
			}
		}
		return n
	}()...)

	c := vgsvg.New(vg.Centimeters(19), vg.Centimeters(10))
	da := plot.MakeDrawArea(c)
	trX, _ := p.Transforms(&da)
	w := ((trX(float64(2*maxLength)) - trX(float64(0))) / vg.Length(2*maxLength)) / 3

	barsFivePrime.Width = w
	barsFivePrime.Offset = -w / 2
	barsThreePrime.Width = w
	barsThreePrime.Offset = w / 2

	p.Draw(da)

	f, err := os.Create(decorate(path, "barchart.svg", filter))
	if err != nil {
		return err
	}
	defer f.Close()
	_, err = c.WriteTo(f)

	return err
}

func boxplot(path string, sets []set) error {
	var (
		fiveEnds  = make([]plotter.Values, len(sets))
		threeEnds = make([]plotter.Values, len(sets))

		err error

		ln int
	)
	for i := range sets {
		fiveEnds[i], err = plotter.CopyValues(&normalised{vals: sets[i].fiveEnd})
		if err != nil {
			return err
		}
		threeEnds[i], err = plotter.CopyValues(&normalised{vals: sets[i].threeEnd})
		if err != nil {
			return err
		}
		if i == 0 {
			ln = fiveEnds[i].Len()
		}
		if fiveEnds[i].Len() != threeEnds[i].Len() || fiveEnds[i].Len() != ln {
			return errors.New("missing values")
		}
	}

	font, err := vg.MakeFont("Helvetica", 10)
	if err != nil {
		return err
	}
	titleFont, err := vg.MakeFont("Helvetica", 12)
	if err != nil {
		return err
	}
	style := plot.TextStyle{Color: color.Gray{0}, Font: font}
	p, err := plot.New()
	if err != nil {
		return err
	}
	p.Title.Text = titles[filter]
	p.Title.TextStyle = plot.TextStyle{Color: color.Gray{0}, Font: titleFont}
	p.X.Label.Text = "Length Offset"
	p.Y.Label.Text = "Relative Frequency"
	p.X.Label.TextStyle = style
	p.Y.Label.TextStyle = style
	p.X.Tick.Label = style
	p.Y.Tick.Label = style
	p.Legend.TextStyle = style

	type boxPair struct{ five, three *plotter.BoxPlot }
	var boxes []boxPair
	for i := 0; i < ln; i++ {
		fiveEnd := make(plotter.Values, len(sets))
		threeEnd := make(plotter.Values, len(sets))
		for j := range sets {
			fiveEnd[j] = fiveEnds[j][i]
			threeEnd[j] = threeEnds[j][i]
		}

		boxFivePrime, err := plotter.NewBoxPlot(1, float64(i), fiveEnd) // A non-zero width is required to prevent the creation failing.
		if err != nil {
			return err
		}
		boxFivePrime.MedianStyle.Width = 0.5
		boxFivePrime.BoxStyle.Width = 0.75
		boxFivePrime.BoxStyle.Color = plotutil.Color(0)

		boxThreePrime, err := plotter.NewBoxPlot(1, float64(i), threeEnd) // A non-zero width is required to prevent the creation failing.
		if err != nil {
			return err
		}
		boxThreePrime.MedianStyle.Width = 0.5
		boxThreePrime.BoxStyle.Width = 0.75
		boxThreePrime.BoxStyle.Color = plotutil.Color(1)

		boxes = append(boxes, boxPair{boxFivePrime, boxThreePrime})

		p.Add(boxFivePrime, boxThreePrime)
	}

	p.Legend.Add("5'-end", &plotter.BarChart{Color: plotutil.Color(0)})
	p.Legend.Add("3'-end", &plotter.BarChart{Color: plotutil.Color(1)})
	p.Legend.Top = true
	p.NominalX(func() []string {
		n := make([]string, ln)
		for i := 0; i < ln; i++ {
			if v := i - maxLength; v%5 == 0 {
				n[i] = fmt.Sprint(v)
			}
		}
		return n
	}()...)
	p.X.Width = 0.5
	p.X.Tick.Width = 0.5
	p.X.Tick.Length = 8
	p.Add(&plotter.Grid{Vertical: plotter.DefaultGridLineStyle})

	c := vgsvg.New(vg.Centimeters(19), vg.Centimeters(10))
	da := plot.MakeDrawArea(c)
	trX, _ := p.Transforms(&da)
	w := ((trX(float64(2*maxLength)) - trX(float64(0))) / vg.Length(2*maxLength)) / 3

	for _, b := range boxes {
		b.five.Width = w
		b.five.Offset = -w / 2
		b.three.Width = w
		b.three.Offset = w / 2
	}

	p.Draw(da)

	f, err := os.Create(decorate(path, "boxplot.svg", filter))
	if err != nil {
		return err
	}
	defer f.Close()
	_, err = c.WriteTo(f)

	return err
}
