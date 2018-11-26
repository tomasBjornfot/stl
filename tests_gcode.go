package main

import (
	"fmt" 
	"io/ioutil"
	"strings"
	"strconv"
	"math"
)
type Points2d struct {
	points [][2]float64
	k	[]float64
	m	[]float64
	d2	[]float64
}

func load2dpoints(path string) [][2]float64 {
	bytes, err := ioutil.ReadFile(path)
	if err != nil {
		fmt.Println("Somthing whent wrong reading the file")
	}
	lines := strings.Split(string(bytes), "\n")
	data := make([][2]float64, len(lines))
	for i := range data {
		s := strings.Split(lines[i], " ")
		data[i][0], _ = strconv.ParseFloat(s[0], 64)
		data[i][1], _ = strconv.ParseFloat(s[1], 64)
	}
	return data
}
func line_square_length(p0 [2]float64, p1 [2]float64) float64 {
	dx := p0[0] - p1[0]
	dy := p0[1] - p1[1]
	return dx*dx + dy*dy
}
func line_k(p0 [2]float64, p1 [2]float64) float64 {
	dx := p1[0] - p0[0]
	dy := p1[1] - p0[1]
	return dy/dx
}
func line_m(p0 [2]float64, p1 [2]float64) float64 {
	k := line_k(p0, p1)
	return p0[1] - k*p0[0]
}
func load_data(path string) *Points2d {
	x := new(Points2d)
	x.points = load2dpoints(path)
	x.d2 = make([]float64, len(x.points)/2)
	x.k = make([]float64, len(x.points)/2)
	x.m = make([]float64, len(x.points)/2)
	for i := range x.d2 {
		x.d2[i] = line_square_length(x.points[2*i], x.points[2*i+1])
		x.k[i] = line_k(x.points[2*i], x.points[2*i+1])
		x.m[i] = line_m(x.points[2*i], x.points[2*i+1])
	}
	return x
}
func print_struct(points *Points2d) {
	for i := range points.d2 {
		fmt.Print("point: ")
		fmt.Printf("%.2f", points.points[2*i])
		fmt.Printf("%.2f", points.points[2*i+1])
		fmt.Print(" d: ")
		fmt.Printf("%.2f", math.Sqrt(points.d2[i]))
		fmt.Print(" k: ")
		fmt.Printf("%.2f", points.k[i])
		fmt.Print(" m: ")
		fmt.Printf("%.2f", points.m[i])
		fmt.Print("\n")
	}
}
func main() {
	
	// populating the struct
	cs := load_data("c:\\tmp\\cs_deck.txt")
	print_struct(cs)
}
