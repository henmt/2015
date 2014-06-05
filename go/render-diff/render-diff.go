// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"code.google.com/p/biogo.graphics/palette/brewer"
	"code.google.com/p/biogo.graphics/rings"
	"code.google.com/p/biogo.rnaseq/norm"
	"code.google.com/p/biogo/feat"
	"code.google.com/p/biogo/feat/genome"
	"code.google.com/p/biogo/feat/genome/mouse/mm10"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/vg"

	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"image/color"
	"math"
	"os"
	"path/filepath"
	"strings"
)

var (
	in     string
	format string

	minLength, maxLength, binLength int

	normalisation int

	minTrace  float64
	maxTrace  float64
	maxCounts float64

	highlight set
	palname   string
)

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
	*s = append(*s, strings.Split(strings.ToLower(value), ",")...)
	if len(*s) == 0 {
		return errors.New("empty set")
	}
	return nil
}

func init() {
	flag.StringVar(&in, "in", "", "json file to be rendered.")
	flag.StringVar(&format, "format", "svg", "specifies the output format of the figure: eps, jpg, jpeg, pdf, png, svg, and tiff.")
	flag.Var(&highlight, "highlight", "comma separated set of chromosome names to highlight.")
	flag.StringVar(&palname, "palette", "Set1", "specify the palette name for highlighting.")
	flag.Float64Var(&minTrace, "tracemin", 0, "set the minimum value for the outer trace if not zero.")
	flag.Float64Var(&maxTrace, "tracemax", 0, "set the maximum value for the outer trace if not zero.")
	flag.Float64Var(&maxCounts, "countmax", 0, "set the maximum value for the inner trace if not zero.")
	flag.IntVar(&normalisation, "normalise", 0, "normalisation strategy: 0 - lib size, 1 - per bin, 2 - per bin/size.")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if in == "" {
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

var index = map[string]int{}

func init() {
	for i, c := range mm10.Chromosomes {
		index[strings.ToLower(c.Chr)] = i
	}
}

func main() {
	rna, err := readJSON(in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	minLength, maxLength, binLength = rna.Min, rna.Max, rna.Bin

	switch normalisation {
	case 0:
		weightFactors[0], weightFactors[1] = float64(rna.Totals[0]), float64(rna.Totals[1])
	case 1:
		weightFactors, err = normaliseByBin(rna)
	case 2:
		weightFactors, err = normaliseByBlock(rna)
	default:
		err = errors.New("illegal normalisation strategy")
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p, err := plot.New()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	mm, lo, hi, err := mouseTracks(rna.Features, highlight, palname, vg.Centimeters(15), rna.Min-rna.Max)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p.Add(mm...)

	p.HideAxes()

	font, err := vg.MakeFont("Helvetica", 14)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	for i, n := range rna.Pair {
		n = filepath.Base(n)
		rna.Pair[i] = n[:strings.Index(n, filepath.Ext(n))]
	}
	p.Title.Text = fmt.Sprintf(
		`%s
%s
%s
min base quality: %v, minimum mapping score: %d
minimum identity: %d%%
length range: [%d,%d]
heat range: [%f,%f]`,
		decorate(rna.Pair[1], format, rna.Filter),
		rna.Classes,
		rna.Pair,
		rna.MinQ, rna.MapQ,
		rna.MinID,
		rna.Min, rna.Max,
		lo, hi)
	p.Title.TextStyle = plot.TextStyle{Color: color.Gray{0}, Font: font}

	for i, class := range rna.Classes {
		rna.Classes[i] = class[strings.LastIndex(class, "/")+1:]
	}

	var classes string
	if len(rna.Classes) > 0 {
		classes = "-" + strings.Join(rna.Classes, ",")
	}
	err = p.Save(vg.Centimeters(19).Inches(), vg.Centimeters(25).Inches(),
		decorate(filepath.Base(rna.Pair[1])+classes, format, rna.Filter),
	)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func decorate(out, format string, filter int) string {
	switch filter {
	case all:
		return fmt.Sprintf("%s-diff.%s", out, format)
	case primary:
		return fmt.Sprintf("%s-U1-diff.%s", out, format)
	case secondary:
		return fmt.Sprintf("%s-A10-diff.%s", out, format)
	default:
		panic("illegal filter")
	}
}

func sum(f []float64) float64 {
	var s float64
	for _, v := range f {
		s += v
	}
	return s
}

func normaliseByBin(rna *Ranged) ([2]float64, error) {
	var data [2][]float64
	for _, f := range rna.Features {
		f := f.(*feature)
		for i := range data {
			data[i] = append(data[i], sum(f.counts[i]))
		}
	}
	var factors [2]float64
	f, err := norm.TMM(data[:], -1, 0.3, 0.05, -1e10, true)
	copy(factors[:], f)
	return factors, err
}

func normaliseByBlock(rna *Ranged) ([2]float64, error) {
	var data [2][]float64
	for _, f := range rna.Features {
		f := f.(*feature)
		for i := range data {
			data[i] = append(data[i], f.counts[i]...)
		}
	}
	var factors [2]float64
	f, err := norm.TMM(data[:], -1, 0.3, 0.05, -1e10, true)
	copy(factors[:], f)
	return factors, err
}

type Ranged struct {
	Pair [2]string

	Bin     int
	Classes []string
	Filter  int

	Min int
	Max int

	MinQ   int
	MinAvQ float64
	MinID  int
	MapQ   int

	Totals   [2]int
	Features []rings.Scorer
}

func readJSON(in string) (rf *Ranged, err error) {
	jsf, err := os.Open(in)
	if err != nil {
		return nil, err
	}
	defer jsf.Close()

	type rangedJSONFeatures struct {
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

	var v rangedJSONFeatures

	err = json.NewDecoder(jsf).Decode(&v)
	if err != nil {
		return nil, err
	}

	rf = &Ranged{
		Pair:    v.Pair,
		Bin:     v.Bin,
		Classes: v.Classes,
		Filter:  v.Filter,
		Min:     v.Min,
		Max:     v.Max,
		MinQ:    v.MinQ,
		MinAvQ:  v.MinAvQ,
		MinID:   v.MinID,
		MapQ:    v.MapQ,
		Totals:  v.Totals,
	}
	rf.Features = make([]rings.Scorer, len(v.Features))
	for i, f := range v.Features {
		rf.Features[i] = f
	}

	return rf, nil
}

var weightFactors [2]float64

type feature struct {
	chr *genome.Chromosome
	start,
	end int

	typ string

	counts   [2][]float64
	supports [2]int
}

func (f *feature) UnmarshalJSON(b []byte) error {
	type jsonFeature struct {
		Chr      string       `json:"chr"`
		Start    int          `json:"start"`
		End      int          `json:"end"`
		Type     string       `json:"type"`
		Counts   [2][]float64 `json:"counts"`
		Supports [2]int       `json:"support"`
	}

	var jf jsonFeature
	err := json.Unmarshal(b, &jf)
	if err != nil {
		return err
	}
	*f = feature{
		chr:      mm10.Chromosomes[index[strings.ToLower(jf.Chr)]],
		start:    jf.Start,
		end:      jf.End,
		typ:      jf.Type,
		counts:   jf.Counts,
		supports: jf.Supports,
	}
	return nil
}

func (f *feature) Start() int { return f.start }
func (f *feature) End() int   { return f.end }
func (f *feature) Len() int   { return f.end - f.start }
func (f *feature) Name() string {
	return fmt.Sprintf("%s[%d,%d)", f.chr.Name(), f.start, f.end)
}
func (f *feature) Description() string    { return f.typ }
func (f *feature) Location() feat.Feature { return f.chr }
func (f *feature) Scores() []float64 {
	scores := make([]float64, len(f.counts[0]))
	for i := range scores {
		scores[i] = f.counts[1][i]/weightFactors[1] - f.counts[0][i]/weightFactors[0]
	}
	return scores
}

type tfs struct{ *feature }

const (
	lomin = 23
	lomax = 28
	himin = 28
	himax = 33
)

func (f tfs) Scores() []float64 {
	var t [2]float64
	scores := f.feature.Scores()
	if lomin-minLength >= 0 && himax-minLength < len(scores) {
		for _, v := range scores[lomin-minLength : lomax-minLength] {
			t[0] += v
		}
		for _, v := range scores[himin-minLength : himax-minLength] {
			t[1] += v
		}
	} else {
		t = [2]float64{math.NaN(), math.NaN()}
	}
	return t[:]
}

type ctfs struct{ *feature }

func (f ctfs) Scores() []float64 {
	factor := float64(binLength) / float64(f.Len())
	return []float64{float64(f.supports[0]) * factor, float64(f.supports[1]) * factor}
}

type symmetricHeat struct{ *rings.Heat }

func (h symmetricHeat) Configure(da plot.DrawArea, cen plot.Point, _ rings.ArcOfer, inner, outer vg.Length, min, max float64) {
	h.Heat.Configure(da, cen, nil, inner, outer, min, max)
	mag := math.Max(math.Abs(h.Min), math.Abs(h.Max))
	if mag == -mag {
		mag++
	}
	h.Min, h.Max = -mag, mag
}

func mouseTracks(scores []rings.Scorer, highlight []string, palname string, diameter vg.Length, lenRange int) (pp []plot.Plotter, lo, hi float64, err error) {
	var p []plot.Plotter

	radius := diameter / 2

	// Relative sizes.
	const (
		gap = 0.005

		label = 117. / 110.

		countsInner = 25. / 110.
		countsOuter = 40. / 110.

		heatInner = 45. / 110.
		heatOuter = 75. / 110.

		traceInner = 80. / 110.
		traceOuter = 95. / 110.

		karyotypeInner = 100. / 110.
		karyotypeOuter = 1.

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
		radius*karyotypeInner, radius*karyotypeOuter, gap,
	)
	if err != nil {
		return nil, 0, 0, err
	}
	mm.LineStyle = sty

	pal, err := brewer.GetPalette(brewer.TypeQualitative, palname, len(highlight))
	if err == nil {
		for i, hn := range highlight {
			for _, c := range mm.Set {
				if hn == strings.ToLower(c.Name()) {
					arc, err := mm.Base.ArcOf(c, nil)
					arc.Theta += rings.Complete * gap / 2
					arc.Phi -= rings.Complete * gap
					if err != nil {
						fmt.Printf("could not find: %s\n", hn)
						break
					}
					col := pal.Colors()[i]
					nc := color.NRGBAModel.Convert(col).(color.NRGBA)
					nc.A = 0x40
					h := rings.NewHighlight(
						nc,
						arc,
						radius*(traceInner-2.5/110.),
						radius*(label+5./110.),
					)
					h.LineStyle = sty
					h.LineStyle.Width /= 4
					p = append(p, h)
					break
				}
			}
		}
	} else if len(highlight) > 0 {
		fmt.Println("no palette")
	}

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
		return nil, 0, 0, fmt.Errorf("bands: %v", err)
	}
	p = append(p, b)
	c, err := rings.NewBlocks(cens, mm, radius*karyotypeInner, radius*karyotypeOuter)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("centromeres: %v", err)
	}
	p = append(p, c)

	font, err := vg.MakeFont("Helvetica", radius*large)
	if err != nil {
		return nil, 0, 0, err
	}
	lb, err := rings.NewLabels(mm, radius*label, rings.NameLabels(mm.Set)...)
	if err != nil {
		return nil, 0, 0, err
	}
	lb.TextStyle = plot.TextStyle{Color: color.Gray16{0}, Font: font}
	p = append(p, lb)

	s, err := rings.NewScores(scores, mm, radius*heatInner, radius*heatOuter,
		symmetricHeat{&rings.Heat{Palette: brewer.Spectral[11].Colors()}},
	)
	if err != nil {
		return nil, 0, 0, err
	}
	p = append(p, s)

	smallFont, err := vg.MakeFont("Helvetica", radius*small)
	if err != nil {
		return nil, 0, 0, err
	}

	traces := make([]rings.Scorer, len(scores))
	for i, s := range scores {
		traces[i] = tfs{s.(*feature)}
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
		return nil, 0, 0, err
	}
	if minTrace != 0 {
		t.Min = minTrace
	}
	if maxTrace != 0 {
		t.Max = maxTrace
	}
	if !math.IsInf(t.Max-t.Min, 0) {
		p = append(p, t)
	}

	counts := make([]rings.Scorer, len(scores))
	for i, s := range scores {
		counts[i] = ctfs{s.(*feature)}
	}
	ct, err := rings.NewScores(counts, mm, radius*countsInner, radius*countsOuter,
		&rings.Trace{
			LineStyles: func() []plot.LineStyle {
				ls := []plot.LineStyle{sty, sty}
				for i, c := range []color.Color{brewer.Set1[4].Colors()[2], brewer.Set1[4].Colors()[3]} {
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
		return nil, 0, 0, err
	}
	if maxCounts != 0 {
		ct.Max = maxCounts
	}
	p = append(p, ct)

	return p, s.Min, s.Max, nil
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
