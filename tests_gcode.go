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
// read/write
func read2dpoints(path string) [][2]float64 {
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
func write2dpoints(path string, points [][2]float64) {
	data_string := " "
	for _, val := range points {
		data_string += fmt.Sprintf("%.2f", val[0])
		data_string += " "
		data_string += fmt.Sprintf("%.2f", val[1])
		data_string += "\n"
	}
	
	bytes := []byte(data_string)
	ioutil.WriteFile(path, bytes, 0666)
}


// for the struct
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
func load_struct_from_file(path string) *Points2d {
	x := new(Points2d)
	x.points = read2dpoints(path)
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
func load_struct(points [][2]float64) *Points2d {
	x := new(Points2d)
	x.points = points
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

func calc_next_point(p_left [2]float64, r2 float64, k float64, m float64) [2]float64 {
	dx := math.Sqrt(r2/(1+k*k))
	p_new := [2]float64 {0.0, 0.0}
	p_new[0] = p_left[0] + dx
	p_new[1] = k*p_new[0] + m
	return p_new
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
	cs := load_struct_from_file("c:\\tmp\\cs_deck.txt")
	//print_struct(cs)
	
	// calculates the new points with a fixed distance
	space := 1.0
	d2_left := 0.0
	p := make([][2]float64, 1000)
	no_p := 0
	line_index := 0
	
	// setting the first value to the same as the first value in cs
	p[0] = cs.points[0] 
	
	for i:=1; i<30; i++ {
		space2 := space*space
		// checks if the next point should be on the "first" line
		d2_left = line_square_length(p[no_p], cs.points[2*line_index+1])
		if space2 < d2_left {
			new_point := calc_next_point(p[no_p], space2, cs.k[line_index], cs.m[line_index])
			p[i] = new_point
			no_p++
			continue
		}
		// checks if the next point should be on the "next" line
		for {
			space2 -= d2_left
			line_index++
			d2_left = cs.d2[line_index]
			if space2 < d2_left {
				new_point := calc_next_point(cs.points[2*line_index], space2, cs.k[line_index], cs.m[line_index])
				p[i] = new_point
				no_p++
				break
			}
		}
	}
	write2dpoints("c:\\tmp\\result.txt", p[:no_p])
}
