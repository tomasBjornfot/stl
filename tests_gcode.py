import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pdb

# This code reads a cross-section generated from an a STL file and adds
# points with equal spaced distance.
# The script assumees that the cross-section points are sorted written 
# as line-pairs (p0 and p1). 

# secure by sorting each line segment

# get line properties
def line_sqr_length(p0, p1):
	dx = p0[0] - p1[0]
	dy = p0[1] - p1[1]
	return dx*dx + dy*dy

def line_k(p0, p1):
	dx = p1[0] - p0[0]
	dy = p1[1] - p0[1]
	k = dy/dx
	return k
	
def line_m(_k, p0):
	return p0[1] - _k*p0[0]

# point properties
def x_dist_to_next_point(r2, k):
	return np.sqrt(r2/(1+k**2))
	
def calc_next(p_left, _r2, _k, _m):
	dx = x_dist_to_next_point(_r2, _k)
	p_new = np.array([0.0, 0.0])
	p_new[0] = p_left[0] + dx
	p_new[1] = _k*p_new[0] + _m
	return p_new

def make_even_spaced(cs, space):
	r2, k, m = [], [], []
	for i in range(len(cs)-1):
		# calculating the s square length of each line segment i cs
		r2.append(line_sqr_length(cs[i], cs[i+1]))
		# calculating the k value for all lines
		k.append(line_k(cs[i], cs[i+1]))
		# calculating the m value for all lines
		m.append(line_m(k[-1], cs[i]))
	# --- generating the new points --- #
	p = []
	# take the first point in cs as a starting point
	p.append(np.array([cs[0, 0], cs[0, 1]]))
	# --- adding points to p that has space2d length == space2 --- #
	line_index = 0
	for i in range(1000):
		#check if new point pass the center line (x=0)
		if p[-1][0] > 0:
			# move the last point to the center line x==0
			p[-1][0] = 0
			p[-1][1] = k[line_index]*0 + m[line_index]
			break
			
		space2 = space**2
		# calculating the square distance from the previous point ...
		# to the end of the line that it belongs to (r2_left)
		r2_left = line_sqr_length(p[-1], cs[line_index+1, :])
		if space2 < r2_left:
			print('i:', i,'first')
			# ... then the next point should be in this line
			# finds the point where the distance is space2 
			# from the last point to the new point
			new_point = calc_next(p[-1], space2, k[line_index], m[line_index])
			p.append(new_point)
		else:
			while True:
				# new space2 when looking at the next line...
				# space2 == space2 - left_r2
				space2 = (np.sqrt(space2) - np.sqrt(r2_left))*(np.sqrt(space2) - np.sqrt(r2_left))
				# checks if there is a new line
				"""
				if line_index + 1 == len(k):
					# no new line
					break
				"""
				# move index to next line
				line_index = line_index + 1

				# move to next line distance. Note that r2_left is the whole
				# line length since previous point was on the previous line
				r2_left = r2[line_index]
				# checks if the point should be in the next line
				if space2 < r2_left:
					print('i:', i,'next')
					print('space2:', space2)
					# ... then the next point should be in this line
					# finds the point where the distance is space2 
					new_point = calc_next(cs[line_index, :], space2, k[line_index], m[line_index])
					p.append(new_point)	
					break
	#pdb.set_trace()
	return p

def unique(data):
	index_remove = []
	for i in range(len(data)-1):
		arr = data[i, :]
		for j in range(i+1, len(data)):
			dx = np.round(arr[0] - data[j,0], 2) 
			if dx == 0:
				index_remove.append(j)
	
	data = np.delete(data, index_remove, 0)			
	return data

def load_cs(row):
	with open('C:\\Go\\data\\_deck_x', 'r') as f:
		lines = f.readlines()
		x = np.array([float(sval) for sval in lines[row].split()])
	with open('C:\\Go\\data\\_deck_z', 'r') as f:
		lines = f.readlines()
		y = np.array([float(sval) for sval in lines[row].split()])
	cs = np.zeros((len(x), 2))
	for i in range(len(x)):
		cs[i, 0] = x[i]
		cs[i, 1] = y[i]
	index = np.argsort(cs[:, 0]) # viktig
	cs = cs[index, :]
	return unique(cs) # viktig
		
# --- MAIN SCRIPT --- #
#cs = np.loadtxt('c:\\tmp\\cs_deck.txt', delimiter=' ')
"""
for i in range(160,165):
	cs = load_cs(i)
	#pdb.set_trace()
	p = make_even_spaced(cs, 4.0)

	plt.figure(figsize=(10,5))
	plt.plot(cs[:, 0], cs[:, 1], '-o', color='C0')
	p = np.array(p)
	plt.plot(p[:, 0], p[:, 1], 'o', color='C1')
	plt.grid()
	plt.title('index: ' + str(i))
	plt.axis('equal')
	plt.show()
"""
# checking the distance between points
"""
d = np.diff(p, axis=0)
distance = np.sqrt( d[:,0]**2 + d[:, 1]**2 )
"""

# --- PLOT EVEN SPACED DATA ---
# reads the y values from _deck_y

y = []
with open('c:\\Go\\data\\_deck_y', 'r') as f:
	lines = f.readlines()
	val = [float(line.split()[0]) for line in lines]
	y.append(val)
y = np.array(y)[0]

# reads the x and z values
x, z = [], []
for i in range(165):
	with open('c:\\Go\\data\\cs_'+str(i)+'.txt', 'r') as f:
		lines = f.readlines()
		x.append(np.array([float(line.split()[0]) for line in lines]))
		z.append(np.array([float(line.split()[1]) for line in lines]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(165):
	y_array = np.zeros(len(x[i])) + y[i]
	ax.scatter(x[i], y_array, z[i])
plt.show()

