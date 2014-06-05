// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"os"
)

func main() {
	fmt.Println(os.Args[1:])
	var (
		minTrace,
		maxTrace float64
		maxCount int
	)
	for i, in := range os.Args[1:] {
		mint, maxt, maxc, err := readJSON(in)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		if i == 0 {
			minTrace = mint
			maxTrace = maxt
			maxCount = maxc
			continue
		}
		minTrace = math.Min(minTrace, mint)
		maxTrace = math.Max(maxTrace, maxt)
		maxCount = max(maxCount, maxc)
	}
	fmt.Printf("Trace minimum: %f\nTrace maximum: %f\nCount maximum: %d\n", minTrace, maxTrace, maxCount)
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

const (
	lomin = 23
	lomax = 28
	himin = 28
	himax = 33
)

func scores(fs []float64, minLength int) (short, long float64) {
	if lomin-minLength >= 0 && himax-minLength < len(fs) {
		for _, v := range fs[lomin-minLength : lomax-minLength] {
			short += v
		}
		for _, v := range fs[himin-minLength : himax-minLength] {
			long += v
		}
		return short, long
	}
	return math.NaN(), math.NaN()
}

func readJSON(in string) (minTrace, maxTrace float64, maxCounts int, err error) {
	jsf, err := os.Open(in)
	if err != nil {
		return math.NaN(), math.NaN(), 0, err
	}
	defer jsf.Close()

	type rangedJSONFeatures struct {
		Min      int               `json:"min"`
		Max      int               `json:"max"`
		Features []json.RawMessage `json:"features"`
	}
	var v rangedJSONFeatures

	err = json.NewDecoder(jsf).Decode(&v)
	if err != nil {
		return math.NaN(), math.NaN(), 0, err
	}
	if len(v.Features) == 0 {
		return math.NaN(), math.NaN(), 0, errors.New("no feature")
	}

	type (
		jsonFeature struct {
			Scores   []float64 `json:"scores"`
			Supports int       `json:"support"`
		}
		jsonFeature2 struct {
			Scores   []float64 `json:"scores"`
			Supports [2]int    `json:"support"`
		}
	)

	var (
		v1 = &jsonFeature{}
		v2 = &jsonFeature2{}

		t interface{}
	)
	if err = json.Unmarshal(v.Features[0], v1); err == nil {
		t = &jsonFeature{}
		maxCounts = v1.Supports
		maxTrace = math.Max(scores(v1.Scores, v.Min))
	} else if err = json.Unmarshal(v.Features[0], v2); err == nil {
		t = &jsonFeature2{}
		maxCounts = max(v2.Supports[0], v2.Supports[1])
		sh, lo := scores(v2.Scores, v.Min)
		minTrace = math.Min(sh, lo)
		maxTrace = math.Max(sh, lo)
	} else {
		return math.NaN(), math.NaN(), 0, err
	}

	for _, f := range v.Features {
		err = json.Unmarshal(f, t)
		if err != nil {
			return math.NaN(), math.NaN(), 0, err
		}
		switch t := t.(type) {
		case *jsonFeature:
			maxCounts = max(maxCounts, t.Supports)
			m := math.Max(scores(t.Scores, v.Min))
			if math.IsNaN(m) {
				break
			}
			maxTrace = math.Max(maxTrace, m)
		case *jsonFeature2:
			maxCounts = max(maxCounts, max(t.Supports[0], t.Supports[1]))
			sh, lo := scores(t.Scores, v.Min)
			min := math.Min(sh, lo)
			max := math.Max(sh, lo)
			if !math.IsNaN(min) {
				minTrace = math.Min(minTrace, min)
			}
			if !math.IsNaN(max) {
				maxTrace = math.Max(maxTrace, max)
			}
		}
	}

	return minTrace, maxTrace, maxCounts, nil
}
