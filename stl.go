package stl

import (
	"os"
	"fmt"
	"math"
	"strings"
	"strconv"
	"encoding/binary"
	"bytes"
)
/*
 * STUCTS
 */
type Mesh struct {
	Triangles [][9]float64
	Normals   [][3]float64
	Profile   [][2]float64
	No_tri    int
	X_min     float64
	Y_min     float64
	Z_min     float64
	X_max     float64
	Y_max     float64
	Z_max     float64
}

/*
 * PRIVATE FUNCTIONS/METHODS
 */
func (mesh *Mesh) calculateMeshProperties() {
	// tar reda på max och min i varje dimension
	mesh.X_min, mesh.Y_min, mesh.Z_min = 100000.0, 100000.0, 100000.0
	mesh.X_max, mesh.Y_max, mesh.Z_max = -100000.0, -100000.0, -100000.0
	points := TrianglesToPoints(*mesh)

	for i := 0; i < 3*mesh.No_tri; i++ {
		// min x
		if points[i][0] < mesh.X_min {
			mesh.X_min = points[i][0]
		}
		// min y
		if points[i][1] < mesh.Y_min {
			mesh.Y_min = points[i][1]
		}
		// min z
		if points[i][2] < mesh.Z_min {
			mesh.Z_min = points[i][2]
		}
		// max x
		if points[i][0] > mesh.X_max {
			mesh.X_max = points[i][0]
		}
		// max y
		if points[i][1] > mesh.Y_max {
			mesh.Y_max = points[i][1]
		}
		// max z
		if points[i][2] > mesh.Z_max {
			mesh.Z_max = points[i][2]
		}
	}
}
func (mesh *Mesh) calculateProfile(radius float64, resolution int) {
	// räknar ut Profilen på brädan i xy planet

	// plockar ut punkterna från trianglarna
	points := TrianglesToPoints(*mesh)

	x := make([]float64, len(points))
	y := make([]float64, len(points))
	for i := 0; i < len(points); i++ {
		x[i] = points[i][1]
		y[i] = points[i][0]
	}

	// tar reda på y värden på Profilen genom att:
	// 0. skapa x värden som cikeln ska ramla ner på
	// 1. iterera över alla x drop värden
	// 2. iterera över alla punkter
	// 3. väljer bara punkter där y > 0
	// 4. tar ut dom punkter som ligger inom [x - radius,x + radius]
	// 5. beräkna höjd på cirkel som krockar med värdet
	// 6. ta ut index på punkten (som cirkeln krockar med)
	// 7. om det finns ett större värde, gäller detta

	// 0
	mesh.calculateMeshProperties()
	drop_x := linspace(mesh.Y_min+0.1-radius, mesh.Y_max-0.1+radius, resolution)
	// 1
	r2 := radius * radius
	index := make([]int, len(drop_x))
	for i := 0; i < len(drop_x); i++ {
		ymax := float64(-1)
		index[i] = -1
		// 2
		for j := 0; j < len(x); j++ {
			// 3
			if y[j] > 0 {
				// 4
				if x[j] > drop_x[i]-radius && x[j] < drop_x[i]+radius {
					// 5
					ymax_new := y[j] + math.Sqrt(r2+(x[j]-drop_x[i])*(x[j]-drop_x[i]))
					if ymax_new > ymax {
						// 6
						ymax = ymax_new
						index[i] = j
					}
				}
			}
		}
	}
	no_index := int(0)
	for i := 0; i < len(index); i++ {
		if index[i] != -1 {
			no_index++
		}
	}	
	mesh.Profile = make([][2]float64, no_index)
	pindex := int(0)
	for i := 0; i < len(index); i++ {
		if index[i] != -1 {
			mesh.Profile[pindex][0] = x[index[i]]
			mesh.Profile[pindex][1] = y[index[i]]
			pindex++
		}
	}
}
func (mesh *Mesh) calculateNormals() {
	// räknar ut normaler för trianglar
	var v0, v1 [3]float64
	for i := 0; i < mesh.No_tri; i++ {
		v0[0] = mesh.Triangles[i][0] - mesh.Triangles[i][6]
		v0[1] = mesh.Triangles[i][1] - mesh.Triangles[i][7]
		v0[2] = mesh.Triangles[i][2] - mesh.Triangles[i][8]

		v1[0] = mesh.Triangles[i][0] - mesh.Triangles[i][3]
		v1[1] = mesh.Triangles[i][1] - mesh.Triangles[i][4]
		v1[2] = mesh.Triangles[i][2] - mesh.Triangles[i][5]

		mesh.Normals[i] = crossProduct(v0, v1)
	}
}
func crossProduct(v0 [3]float64, v1 [3]float64) [3]float64 {
	var x [3]float64
	x[0] = v0[1]*v1[2] - v0[2]*v1[1]
	x[1] = v0[2]*v1[0] - v0[0]*v1[2]
	x[2] = v0[0]*v1[1] - v0[1]*v1[0]
	length := math.Sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
	x[0] = x[0] / length
	x[1] = x[1] / length
	x[2] = x[2] / length
	return x
}
func linspace(min float64, max float64, no_segments int) []float64 {
	numbers := make([]float64, no_segments+1)
	numbers[0] = min
	segment := (max - min) / float64(no_segments)
	for i := 1; i < len(numbers); i++ {
		numbers[i] = numbers[i-1] + segment
	}
	return numbers
}
func TrianglesToPoints(mesh Mesh) [][3]float64 {
	// tar ut trianglarna från mesh och gör en x,3 matris av punkterna
	points := make([][3]float64, 3*mesh.No_tri)
	// plockar ut punkterna från trianglarna
	for i := 0; i < mesh.No_tri; i++ {
		points[3*i+0][0] = mesh.Triangles[i][0]
		points[3*i+1][0] = mesh.Triangles[i][3]
		points[3*i+2][0] = mesh.Triangles[i][6]

		points[3*i+0][1] = mesh.Triangles[i][1]
		points[3*i+1][1] = mesh.Triangles[i][4]
		points[3*i+2][1] = mesh.Triangles[i][7]

		points[3*i+0][2] = mesh.Triangles[i][2]
		points[3*i+1][2] = mesh.Triangles[i][5]
		points[3*i+2][2] = mesh.Triangles[i][8]
	}
	return points
}

/*
 * PUBLIC FUNCTIONS/METHODS
 */
func (mesh *Mesh) Read(path string, filetype int) {
	// läser in en STL fil till mesh
	file, err := os.Open(path)
	if err != nil {
		fmt.Printf("ReadFromFile: Something went wrong!!!")
	}
	defer file.Close()
	finfo, _ := file.Stat()

	bytes := make([]byte, finfo.Size())
	file.Read(bytes)

	var row int
	s_array := strings.Split(string(bytes), "\n")
	if filetype == 1 {
		mesh.No_tri = (len(s_array)-2)/7 - 1
		mesh.Triangles = make([][9]float64, mesh.No_tri)
		mesh.Normals = make([][3]float64, mesh.No_tri)
		for i := 1; i < mesh.No_tri+1; i++ {
			s_tri := s_array[7*i+1 : 7*(i+1)+1]
			n_string := strings.Split(s_tri[0], " ")

			// läser in normalerna till triangeln
			for j := 0; j < 3; j++ {
				mesh.Normals[row][j], _ = strconv.ParseFloat(n_string[j+2], 32)
			}
			// läser in hörnpunkterna till triangeln
			t1_string := strings.Split(s_tri[2], " ")
			t2_string := strings.Split(s_tri[4], " ")
			t3_string := strings.Split(s_tri[3], " ")
			for j := 0; j < 3; j++ {
				mesh.Triangles[row][j+0], _ = strconv.ParseFloat(t1_string[j+3], 32)
				mesh.Triangles[row][j+3], _ = strconv.ParseFloat(t2_string[j+3], 32)
				mesh.Triangles[row][j+6], _ = strconv.ParseFloat(t3_string[j+3], 32)
			}
			row++
		}
	}
	if filetype == 0 {
		mesh.No_tri = (len(bytes) - 84) / 50
		mesh.Triangles = make([][9]float64, mesh.No_tri)
		mesh.Normals = make([][3]float64, mesh.No_tri)
		bits := binary.LittleEndian.Uint32(bytes[0:4])
		var start_index int = 80 // hoppar över header
		start_index += 4         //  hoppar över uint32
		for i := 0; i < mesh.No_tri; i++ {
			// normaler
			for j := 0; j < 3; j++ {
				bits = binary.LittleEndian.Uint32(bytes[(start_index):(start_index + 4)])
				mesh.Normals[i][j] = float64(math.Float32frombits(bits))
				start_index += 4
			}
			// trianglar
			for j := 0; j < 9; j++ {
				bits = binary.LittleEndian.Uint32(bytes[(start_index):(start_index + 4)])
				mesh.Triangles[i][j] = float64(math.Float32frombits(bits))
				start_index += 4
			}
			// hoppar över uint16
			start_index += 2
		}
	}
	mesh.calculateProfile(50.0, 100)
}
func (mesh *Mesh) Write(path string) {
	// skriver en mesh till binär STL fil
	file, err := os.Create(path)
	if err != nil {
		fmt.Println("WriteToFile: Something went wrong!!!")
	}
	defer file.Close()

	header := make([]byte, 80)
	file.Write(header)
	buf := new(bytes.Buffer)
	binary.Write(buf, binary.LittleEndian, int32(mesh.No_tri))
	for i := 0; i < mesh.No_tri; i++ {
		for j := 0; j < 3; j++ {
			binary.Write(buf, binary.LittleEndian, float32(mesh.Normals[i][j]))
		}

		for j := 0; j < 9; j++ {
			binary.Write(buf, binary.LittleEndian, float32(mesh.Triangles[i][j]))
		}

		binary.Write(buf, binary.LittleEndian, int16(i))
	}
	buf.WriteTo(file)
}
func (mesh *Mesh) MoveToCenter() {
	// flyttar en mesh till centrum
	var x_sum, y_sum, z_sum float64 = 0, 0, 0
	var vector [3]float64
	for i := 0; i < mesh.No_tri; i++ {
		x_sum += mesh.Triangles[i][0]
		x_sum += mesh.Triangles[i][3]
		x_sum += mesh.Triangles[i][6]

		y_sum += mesh.Triangles[i][1]
		y_sum += mesh.Triangles[i][4]
		y_sum += mesh.Triangles[i][7]

		z_sum += mesh.Triangles[i][2]
		z_sum += mesh.Triangles[i][5]
		z_sum += mesh.Triangles[i][8]
	}
	vector[0] = -x_sum / (3 * float64(mesh.No_tri))
	vector[1] = -y_sum / (3 * float64(mesh.No_tri))
	vector[2] = -z_sum / (3 * float64(mesh.No_tri))
	fmt.Println("translation vector: ", vector)
	mesh.Translate(vector)
}
func (mesh *Mesh) MoveToCenter2() {
	// flyttar till centrum map högsta och minsta värde
	var X_min, Y_min, Z_min float64 = 1000.0, 1000.0, 1000.0
	var X_max, Y_max, Z_max float64 = -1000.0, -1000.0, -1000.0
	var vector [3]float64
	for i := 0; i < mesh.No_tri; i++ {
		for j := 0; j < 7; j = j + 3 {
			if mesh.Triangles[i][j] < X_min {
				X_min = mesh.Triangles[i][j]
			}
			if mesh.Triangles[i][j] > X_max {
				X_max = mesh.Triangles[i][j]
			}
		}
		for j := 1; j < 8; j = j + 3 {
			if mesh.Triangles[i][j] < Y_min {
				Y_min = mesh.Triangles[i][j]
			}
			if mesh.Triangles[i][j] > Y_max {
				Y_max = mesh.Triangles[i][j]
			}
		}
		for j := 2; j < 9; j = j + 3 {
			if mesh.Triangles[i][j] < Z_min {
				Z_min = mesh.Triangles[i][j]
			}
			if mesh.Triangles[i][j] > Z_max {
				Z_max = mesh.Triangles[i][j]
			}
		}
	}
	vector[0] = -(X_max + X_min) / 2.0
	vector[1] = -(Y_max + Y_min) / 2.0
	vector[2] = -(Z_max + Z_min) / 2.0
	fmt.Println("translation vector: ", vector)
	mesh.Translate(vector)
}
func (mesh *Mesh) Translate(vector [3]float64) {
	// translaterar en mesh i rymden
	for i := 0; i < mesh.No_tri; i++ {
		for j := 0; j < 9; j = j + 3 {
			mesh.Triangles[i][j] = mesh.Triangles[i][j] + vector[0]
		}
		for j := 1; j < 9; j = j + 3 {
			mesh.Triangles[i][j] = mesh.Triangles[i][j] + vector[1]
		}
		for j := 2; j < 9; j = j + 3 {
			mesh.Triangles[i][j] = mesh.Triangles[i][j] + vector[2]
		}
	}
}
func (mesh *Mesh) Rotate(axis string, angle_degrees float64) {
	// roterar en mesh i rymden
	ar := angle_degrees * math.Pi / 180.0
	for i := 0; i < mesh.No_tri; i++ {
		for j := 0; j < 9; j = j + 3 {
			point := mesh.Triangles[i][j:(j + 3)]
			new_point := make([]float64, 3)
			if axis == "x" {
				new_point[0] = point[0]
				new_point[1] = point[1]*math.Cos(ar) - point[2]*math.Sin(ar)
				new_point[2] = point[1]*math.Sin(ar) + point[2]*math.Cos(ar)
			}
			if axis == "y" {
				new_point[0] = point[0]*math.Cos(ar) + point[2]*math.Sin(ar)
				new_point[1] = point[1]
				new_point[2] = -point[0]*math.Sin(ar) + point[2]*math.Cos(ar)
			}
			if axis == "z" {
				new_point[0] = point[0]*math.Cos(ar) - point[1]*math.Sin(ar)
				new_point[1] = point[0]*math.Sin(ar) + point[1]*math.Cos(ar)
				new_point[2] = point[2]
			}
			mesh.Triangles[i][j+0] = new_point[0]
			mesh.Triangles[i][j+1] = new_point[1]
			mesh.Triangles[i][j+2] = new_point[2]
		}
	}
	mesh.calculateNormals()
}
func (mesh *Mesh) AlignMesh(cadtype string) {
	// gör en alignment beroende på vilkan cad som används
	if cadtype == "boardcad" {
		mesh.MoveToCenter()
		mesh.Rotate("x", 90)
		mesh.Rotate("z", 90)
	}
	mesh.calculateNormals()
}
func (mesh *Mesh) AlignMeshX() {
	// roterar brädan runt x vectorn tills man hittar ett minimum Z_max - Z_min
	mesh.calculateMeshProperties()
	rmesh := mesh
	z_range := rmesh.Z_max - rmesh.Z_min
	for i := 0; i < 50; i++ {
		rmesh.Rotate("x", -0.1)
		rmesh.calculateMeshProperties()
		if rmesh.Z_max-rmesh.Z_min < z_range {
			z_range = rmesh.Z_max - rmesh.Z_min
		} else {
			//fmt.Printf("Alignment x rotation: %0.2f degrees\n", 0.1*float64(i))
			break
		}
	}
	mesh = rmesh
}
func (mesh *Mesh) Split() (*Mesh, *Mesh) {
	// delar upp deck och bottom på brädan
	// flytta funktionen
	No_tri_deck := int(0)
	No_tri_bottom := int(0)
	for i := 0; i < mesh.No_tri; i++ {
		if mesh.Normals[i][2] < 0 {
			No_tri_deck++
		}
		if mesh.Normals[i][2] >= 0 {
			No_tri_bottom++
		}
	}
	deck := new(Mesh)
	deck.Triangles = make([][9]float64, No_tri_deck)
	deck.Normals = make([][3]float64, No_tri_deck)
	deck.No_tri = No_tri_deck

	bottom := new(Mesh)
	bottom.Triangles = make([][9]float64, No_tri_bottom)
	bottom.Normals = make([][3]float64, No_tri_bottom)
	bottom.No_tri = No_tri_bottom

	i_deck := int(0)
	i_bottom := int(0)
	for i := 0; i < mesh.No_tri; i++ {
		if mesh.Normals[i][2] < 0 {
			deck.Triangles[i_deck] = mesh.Triangles[i]
			deck.Normals[i_deck] = mesh.Normals[i]
			i_deck++
		}
		if mesh.Normals[i][2] >= 0 {
			bottom.Triangles[i_bottom] = mesh.Triangles[i]
			bottom.Normals[i_bottom] = mesh.Normals[i]
			i_bottom++
		}
	}
	deck.calculateProfile(50.0, 100)
	bottom.calculateProfile(50.0, 100)
	return deck, bottom
}

/*
 * EXTRAS
 */
func (mesh *Mesh) WritePointsToFile(path string) {
	// skriver alla triangelhörn på fil
	file, err := os.Create(path)
	if err != nil {
		fmt.Println("WritePointsToFile: Something went wrong!!!")
	}
	for i := 0; i < mesh.No_tri; i++ {
		s_x := strconv.FormatFloat(mesh.Triangles[i][0], 'f', 2, 64)
		s_y := strconv.FormatFloat(mesh.Triangles[i][1], 'f', 2, 64)
		s_z := strconv.FormatFloat(mesh.Triangles[i][2], 'f', 2, 64)
		s := s_x + " " + s_y + " " + s_z + "\n"
		file.WriteString(s)

		s_x = strconv.FormatFloat(mesh.Triangles[i][3], 'f', 2, 64)
		s_y = strconv.FormatFloat(mesh.Triangles[i][4], 'f', 2, 64)
		s_z = strconv.FormatFloat(mesh.Triangles[i][5], 'f', 2, 64)
		s = s_x + " " + s_y + " " + s_z + "\n"
		file.WriteString(s)

		s_x = strconv.FormatFloat(mesh.Triangles[i][6], 'f', 2, 64)
		s_y = strconv.FormatFloat(mesh.Triangles[i][7], 'f', 2, 64)
		s_z = strconv.FormatFloat(mesh.Triangles[i][8], 'f', 2, 64)
		s = s_x + " " + s_y + " " + s_z + "\n"
		file.WriteString(s)
	}
}
func (mesh *Mesh) WriteNormalsToFile(path string) {
	// skriver alla triangelnormaler på fil
	file, err := os.Create(path)
	if err != nil {
		fmt.Println("WriteNormalsToFile: Something went wrong!!!")
	}
	for i := 0; i < mesh.No_tri; i++ {
		s_x := strconv.FormatFloat(mesh.Normals[i][0], 'f', 2, 64)
		s_y := strconv.FormatFloat(mesh.Normals[i][1], 'f', 2, 64)
		s_z := strconv.FormatFloat(mesh.Normals[i][2], 'f', 2, 64)
		s := s_x + " " + s_y + " " + s_z + "\n"
		file.WriteString(s)
	}
}
func (mesh *Mesh) WriteProfileToFile(path string) {
	// Skriver Profilen till fil
	file, err := os.Create(path)
	if err != nil {
		fmt.Println("WriteProfileToFile: Something went wrong!!!")
	}
	for i := 0; i < len(mesh.Profile); i++ {
		s_x := strconv.FormatFloat(mesh.Profile[i][0], 'f', 2, 64)
		s_y := strconv.FormatFloat(mesh.Profile[i][1], 'f', 2, 64)
		s := s_x + " " + s_y + "\n"
		file.WriteString(s)
	}
}
func (mesh *Mesh) WriteMeshProperties() {
	// skriver alla mesh properties på terminal
	mesh.calculateMeshProperties()
	fmt.Println("Mesh properties:")
	fmt.Printf("width x: %.2f mm\n", mesh.X_max-mesh.X_min)
	fmt.Printf("width y: %.2f mm\n", mesh.Y_max-mesh.Y_min)
	fmt.Printf("width z: %.2f mm\n", mesh.Z_max-mesh.Z_min)
}
