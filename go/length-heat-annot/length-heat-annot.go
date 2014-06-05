// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"code.google.com/p/biogo.boom"
	"code.google.com/p/biogo.graphics/palette"
	"code.google.com/p/biogo.graphics/palette/brewer"
	"code.google.com/p/biogo.graphics/rings"
	"code.google.com/p/biogo.store/interval"
	"code.google.com/p/biogo/feat"
	"code.google.com/p/biogo/feat/genome"
	"code.google.com/p/biogo/feat/genome/mouse/mm10"
	"code.google.com/p/biogo/io/featio"
	"code.google.com/p/biogo/io/featio/gff"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/vg"

	"errors"
	"flag"
	"fmt"
	"image/color"
	"io"
	"math"
	"os"
	"strings"
	"unsafe"
)

var (
	in, annot, out string
	classes        set

	filter int

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
	return fmt.Sprint(*s)
}
func (s *set) Set(value string) error {
	*s = append(*s, strings.Split(value, ",")...)
	if len(*s) == 0 {
		return errors.New("set: empty set")
	}
	return nil
}

func init() {
	flag.StringVar(&in, "in", "", "file name of a BAM file to be processed.")
	flag.StringVar(&annot, "annot", "", "file name of a GFF file containing annotations.")
	flag.Var(&classes, "class", "comma separated set of annotation classes to analyse.")
	flag.StringVar(&out, "out", "", "outfile name.")
	flag.IntVar(&minLength, "min", 20, "minimum length read considered.")
	flag.IntVar(&maxLength, "max", 35, "maximum length read considered.")
	flag.IntVar(&minId, "minid", 90, "minimum percentage identity for non-clipped bases.")
	flag.IntVar(&minQ, "minQ", 20, "minimum per-base sequence quality.")
	flag.IntVar(&filter, "f", 0, "filter on piwi type 0: no filter, 1: primary, 2: secondary.")
	flag.Float64Var(&minAvQ, "minAvQ", 30, "minimum average per-base sequence quality.")
	flag.IntVar(&mapQ, "mapQ", 0, "minimum mapping quality [0, 255).")
	flag.IntVar(&binLength, "bin", 1e7, "bin length.")
	flag.StringVar(&format, "format", "svg", "specifies the output format of the example: eps, jpg, jpeg, pdf, png, svg, and tiff.")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	mapQb = byte(mapQ)
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if in == "" || out == "" || annot == "" || len(classes) == 0 || mapQ < 0 || mapQ > 254 {
		flag.Usage()
		os.Exit(1)
	}
	for _, s := range []string{"eps", "jpg", "jpeg", "pdf", "png", "svg", "tiff"} {
		if format == s {
			return
		}
	}
	flag.Usage()
	os.Exit(1)
}

func qualOk(r *boom.Record, minId, minQ int, minAvQ float64) (ok bool) {
	var (
		off, l int
		match  int
		mQ     int
		edit   int
		allOk  = true
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
				allOk = allOk && int(q) >= minQ
				if t == boom.CigarMatch || t == boom.CigarEqual || t == boom.CigarSoftClipped {
					match++
					mQ += int(q)
				}
			}
		}
	}
	match -= edit

	ok = allOk && match*100 >= minId*l && mQ >= int(minAvQ*float64(l))

	return
}

type Location struct {
	Rid int
	Pos int
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

func filterFeats(annot string, classes []string, bam string) ([]interval.IntTree, error) {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		return nil, err
	}
	defer bf.Close()
	names := bf.RefNames()
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

func (f intBam) Range() interval.IntRange { return interval.IntRange{f.Start(), f.End()} }
func (f intBam) Overlap(b interval.IntRange) bool {
	return f.Start()+len(f.Seq()) > b.Start && f.Start() < b.End
}
func (f intBam) ID() uintptr { return *(*uintptr)(unsafe.Pointer(f.Record)) }

func rnaFeats(in string, classFilt []interval.IntTree, minLength, maxLength, minId, minQ int, mapQb byte, minAvQ float64) ([]rings.Scorer, error) {
	bf, err := boom.OpenBAM(in)
	if err != nil {
		return nil, err
	}
	defer bf.Close()

	smap := make(map[Location]map[int]int)
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
			if qualOk(r, minId, minQ, minAvQ) {
				if len(classFilt[r.RefID()].Get(intBam{r})) == 0 {
					continue loop
				}
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
					var sc map[int]int
					loc := Location{Rid: r.RefID(), Pos: r.Start() / binLength}
					sc, ok := smap[loc]
					if !ok {
						sc = make(map[int]int)
					}
					sc[l]++
					smap[loc] = sc
				}
			}
		}
	}
	if err != nil {
		return nil, err
	}

	var scoreFeats []rings.Scorer
	names := bf.RefNames()
	for k, v := range smap {
		c := mm10.Chromosomes[index[strings.ToLower(names[k.Rid])]]
		f := &fs{
			start:    k.Pos * binLength,
			end:      min((k.Pos+1)*binLength, c.Len()),
			location: c,
			scores:   make([]float64, maxLength-minLength+1),
		}
		for l, sv := range v {
			f.scores[l-minLength] = float64(sv) / float64(f.Len())
		}
		scoreFeats = append(scoreFeats, f)
	}

	// Fill any holes due to absence of reads.
	for cid, name := range names {
		c := mm10.Chromosomes[index[strings.ToLower(name)]]
		for p := c.Start(); p < c.End()/binLength; p++ {
			if _, ok := smap[Location{Rid: cid, Pos: p}]; !ok {
				scoreFeats = append(scoreFeats, &fs{
					start:    p * binLength,
					end:      min((p+1)*binLength, c.Len()),
					location: c,
					scores:   make([]float64, maxLength-minLength+1),
				})
			}
		}
	}

	return scoreFeats, nil
}

func main() {
	p, err := plot.New()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	classFilt, err := filterFeats(annot, classes, in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	rna, err := rnaFeats(in, classFilt, minLength, maxLength, minId, minQ, mapQb, minAvQ)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	mm, err := mouseTracks(rna, vg.Centimeters(15), maxLength-minLength)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p.Add(mm...)
	p.HideAxes()

	var name string
	switch filter {
	case all:
		name = fmt.Sprintf("%s.%s", out, format)
	case primary:
		name = fmt.Sprintf("%s-U1.%s", out, format)
	case secondary:
		name = fmt.Sprintf("%s-A10.%s", out, format)
	}
	err = p.Save(vg.Centimeters(19).Inches(), vg.Centimeters(25).Inches(), name)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

type fs struct {
	start, end int
	name       string
	location   feat.Feature
	scores     []float64
}

func (f *fs) Start() int             { return f.start }
func (f *fs) End() int               { return f.end }
func (f *fs) Len() int               { return f.end - f.start }
func (f *fs) Name() string           { return f.name }
func (f *fs) Description() string    { return "alignment bin" }
func (f *fs) Location() feat.Feature { return f.location }
func (f *fs) Scores() []float64      { return f.scores }

type tfs struct{ *fs }

const (
	lomin = 23
	lomax = 28
	himin = 28
	himax = 33
)

func (f *tfs) Scores() []float64 {
	var t [2]float64
	if lomin-minLength >= 0 && himax-minLength < len(f.scores) {
		for _, v := range f.scores[lomin-minLength : lomax-minLength] {
			t[0] += v
		}
		for _, v := range f.scores[himin-minLength : himax-minLength] {
			t[1] += v
		}
	} else {
		t = [2]float64{math.NaN(), math.NaN()}
	}
	return t[:]
}

func mouseTracks(scores []rings.Scorer, diameter vg.Length, lenRange int) ([]plot.Plotter, error) {
	var p []plot.Plotter

	radius := diameter / 2

	// Relative sizes.
	const (
		label = 117. / 110.

		karyotypeInner = 100. / 110.
		karyotypeOuter = 1.

		heatInner = 30. / 110.
		heatOuter = 75. / 110.

		traceInner = 80. / 110.
		traceOuter = 95. / 110.

		large = 7. / 110.
		small = 2. / 110.
	)

	sty := plotter.DefaultLineStyle
	sty.Width /= 2

	chr := make([]feat.Feature, len(mm10.Chromosomes))
	for i, c := range mm10.Chromosomes {
		chr[i] = c
	}
	mm, err := rings.NewGappedBlocks(
		chr,
		rings.Arc{rings.Complete / 4 * rings.CounterClockwise, rings.Complete * rings.Clockwise},
		radius*karyotypeInner, radius*karyotypeOuter, 0.005,
	)
	if err != nil {
		return nil, err
	}
	mm.LineStyle = sty
	p = append(p, mm)

	bands := make([]feat.Feature, len(mm10.Bands))
	cens := make([]feat.Feature, 0, len(mm10.Chromosomes))
	for i, b := range mm10.Bands {
		bands[i] = colorBand{b}
		s := b.Start()
		// This condition depends on p -> q sort order in the $karyotype.Bands variable.
		// All standard genome packages follow this, though here the test is more general than
		// actually required since mm is telocentric.
		if b.Band[0] == 'q' && (s == 0 || mm10.Bands[i-1].Band[0] == 'p') {
			cens = append(cens, colorBand{&genome.Band{Band: "cen", Desc: "Band", StartPos: s, EndPos: s, Giemsa: "acen", Chr: b.Location()}})
		}
	}
	b, err := rings.NewBlocks(bands, mm, radius*karyotypeInner, radius*karyotypeOuter)
	if err != nil {
		return nil, fmt.Errorf("bands: %v", err)
	}
	p = append(p, b)
	c, err := rings.NewBlocks(cens, mm, radius*karyotypeInner, radius*karyotypeOuter)
	if err != nil {
		return nil, fmt.Errorf("centromeres: %v", err)
	}
	p = append(p, c)

	font, err := vg.MakeFont("Helvetica", radius*large)
	if err != nil {
		return nil, err
	}
	lb, err := rings.NewLabels(mm, radius*label, rings.NameLabels(mm.Set)...)
	if err != nil {
		return nil, err
	}
	lb.TextStyle = plot.TextStyle{Color: color.Gray16{0}, Font: font}
	p = append(p, lb)

	s, err := rings.NewScores(scores, mm, radius*heatInner, radius*heatOuter,
		&rings.Heat{Palette: palette.Heat(10, 1).Colors()},
	)
	if err != nil {
		return nil, err
	}
	p = append(p, s)

	traces := make([]rings.Scorer, len(scores))
	for i, s := range scores {
		traces[i] = &tfs{s.(*fs)}
	}

	smallFont, err := vg.MakeFont("Helvetica", radius*small)
	if err != nil {
		return nil, err
	}
	t, err := rings.NewScores(traces, mm, radius*traceInner, radius*traceOuter,
		&rings.Trace{
			LineStyles: func() []plot.LineStyle {
				ls := []plot.LineStyle{sty, sty}
				for i, c := range brewer.Set1[3].Colors()[:len(ls)] {
					nc := color.NRGBAModel.Convert(c).(color.NRGBA)
					nc.A = 0x80
					ls[i].Color = nc
				}
				return ls
			}(),
			Join: true,
			Axis: &rings.Axis{
				Angle:     rings.Complete / 4,
				Grid:      plotter.DefaultGridLineStyle,
				LineStyle: sty,
				Tick: rings.TickConfig{
					Marker:    plot.DefaultTicks,
					LineStyle: sty,
					Length:    2,
					Label:     plot.TextStyle{Color: color.Gray16{0}, Font: smallFont},
				},
			},
		},
	)
	if err != nil {
		return nil, err
	}
	p = append(p, t)

	return p, nil
}

type colorBand struct {
	*genome.Band
}

func (b colorBand) FillColor() color.Color {
	switch b.Giemsa {
	case "acen":
		return color.RGBA{R: 0xff, A: 0xff}
	case "gneg":
		return color.Gray{0xff}
	case "gpos25":
		return color.Gray{3 * math.MaxUint8 / 4}
	case "gpos33":
		return color.Gray{2 * math.MaxUint8 / 3}
	case "gpos50":
		return color.Gray{math.MaxUint8 / 2}
	case "gpos66":
		return color.Gray{math.MaxUint8 / 3}
	case "gpos75":
		return color.Gray{math.MaxUint8 / 4}
	case "gpos100":
		return color.Gray{0x0}
	default:
		panic("unexpected giemsa value")
	}
}

func (b colorBand) LineStyle() plot.LineStyle {
	switch b.Giemsa {
	case "acen":
		return plot.LineStyle{Color: color.RGBA{R: 0xff, A: 0xff}, Width: 1}
	case "gneg", "gpos25", "gpos33", "gpos50", "gpos66", "gpos75", "gpos100":
		return plot.LineStyle{}
	default:
		panic("unexpected giemsa value")
	}
}
