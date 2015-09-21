package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

func mustAtoi(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

type vector []int

func (v *vector) inc(i int) {
	switch {
	case i < len(*v):
		(*v)[i]++
	case i == len(*v):
		*v = append(*v, 1)
	default:
		t := make([]int, i+1)
		copy(t, *v)
		t[i] = 1
		*v = t
	}
}

func main() {
	vm := make(map[string]*vector)
	fm := make(map[string]string)

	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		if f.Feature != "repeat" {
			continue
		}
		att := f.FeatAttributes.Get("repeat")
		if att == "" {
			continue
		}
		fields := strings.Fields(att)

		s := mustAtoi(fields[2])
		e := mustAtoi(fields[3])
		v, ok := vm[fields[0]]
		if !ok {
			v = &vector{}
			vm[fields[0]] = v
			fm[fields[0]] = fields[1]
		}
		for i := e; i >= s; i-- {
			v.inc(i)
		}
	}

	for typ, vec := range vm {
		for pos, val := range *vec {
			if val != 0 {
				fmt.Printf("%s\t%s\t%d\t%d\n", fm[typ], typ, pos, val)
			}
		}
	}
}
