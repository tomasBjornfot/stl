package main

import (
	"fmt" 
	"io/ioutil"
	"strings"
	"strconv"
	"math"
)

// --- useful struct ---
type Points2d struct {
	points [][2]float64
	k	[]float64
	m	[]float64
	d2	[]float64
}
// --- read/write (for testing) ---
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

func read_matrix_row(path string, row int) []float64 {
	bytes, err := ioutil.ReadFile(path)
	if err != nil {
		fmt.Println("Somthing whent wrong reading the file")
	}
	lines := strings.Split(string(bytes), "\n")
	line := strings.Split(lines[row], " ")
	data := make([]float64, len(line))
	for i, val := range line {
		data[i], _ = strconv.ParseFloat(val, 64)
	}
	return data
}

func write2dpoints(path string, points [][2]float64) {
	data_string := ""
	for _, val := range points {
		data_string += fmt.Sprintf("%.2f", val[0])
		data_string += " "
		data_string += fmt.Sprintf("%.2f", val[1])
		data_string += "\n"
	}
	
	bytes := []byte(data_string)
	ioutil.WriteFile(path, bytes, 0666)
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
// --- printing (for testing) ---
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

// --- 2D line calculations ---
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

func load_struct(x []float64, y []float64) *Points2d {
	p2d := new(Points2d)
	p := make([][2]float64, len(x))
	for i := range x {
		p[i][0] = x[i]
		p[i][1] = y[i]
	}
	p2d.points = p
	p2d.d2 = make([]float64, len(p2d.points)/2)
	p2d.k = make([]float64, len(p2d.points)/2)
	p2d.m = make([]float64, len(p2d.points)/2)
	for i := range p2d.d2 {
		p2d.d2[i] = line_square_length(p2d.points[2*i], p2d.points[2*i+1])
		p2d.k[i] = line_k(p2d.points[2*i], p2d.points[2*i+1])
		p2d.m[i] = line_m(p2d.points[2*i], p2d.points[2*i+1])
	}
	return p2d
}

// function to make the even spaced line
func even_spaced_cs(cs *Points2d, space float64) [][2]float64 {
	
	calc_next_point := func(p_left [2]float64, r2 float64, k float64, m float64) [2]float64 {
		dx := math.Sqrt(r2/(1+k*k))
		p_new := [2]float64 {0.0, 0.0}
		p_new[0] = p_left[0] + dx
		p_new[1] = k*p_new[0] + m
		return p_new
	}
	
	// calculates the new points with a fixed distance
	d2_left := 0.0
	p := make([][2]float64, 1000)
	no_p := 0
	line_index := 0
	
	// setting the first value to the same as the first value in cs
	p[0] = cs.points[0] 
	
	for i:=1; i<1000; i++ {
		// checks if new point pass the center line (x=0)
		if p[no_p][0] > 0 {
			// move the point to the center line
			p[no_p][0] = 0
			p[no_p][1] = cs.k[line_index]*p[no_p][0] + cs.m[line_index]
			break
		}
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
			space2 = (math.Sqrt(space2) - math.Sqrt(d2_left))*(math.Sqrt(space2) - math.Sqrt(d2_left))
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
	return p[:no_p+1]
}

func main() {
	// reading a cross-section from a mesh
	// x and z is csMesh.x[row] and csMesh.z[row]
	x := read_matrix_row("c:\\Go\\data\\_deck_x", 1)
	z := read_matrix_row("c:\\Go\\data\\_deck_z", 1)
	// generates a Points2d struct 
	cs := load_struct(x, z)
	// generates an even spaced cross-section
	p := even_spaced_cs(cs, 1.0)
	
	write2dpoints("c:\\Go\\data\\cs_row.txt", cs.points)
	write2dpoints("c:\\Go\\data\\result.txt", p)
}
