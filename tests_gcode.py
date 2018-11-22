import numpy as np
import matplotlib.pyplot as plt
import pdb

# This code reads a cross-section generated from an a STL file and adds
# points with equal spaced distance.
# The script assumees that the cross-section points are sorted written 
# as line-pairs (p0 and p1). 

# get line properties
def line_sqr_length(p0, p1):
	dx = p0[0] - p1[0]
	dz = p0[1] - p1[1]
	return dx*dx + dz*dz

def line_k(p0, p1):
	dx = p1[0] - p0[0]
	dy = p1[1] - p0[1]
	return dy/dx
	
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
	#pdb.set_trace()
	return p_new

# --- MAIN SCRIPT --- #
# --- indata --- #
cs = np.loadtxt('c:\\tmp\\cs_deck.txt', delimiter=' ')
space = 1.0
r2, k, m = [], [], []
for i in range(0,len(cs),2):
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
for i in range(30):
	space2 = space**2
	# calculating the square distance from the previous point ...
	# to the end of the line that it belongs to (r2_left)
	r2_left = line_sqr_length(p[-1], cs[2*line_index+1, :])
	if space2 < r2_left:
		# ... then the next point should be in this line
		# finds the point where the distance is space2 
		# from the last point to the new point
		new_point = calc_next(p[-1], space2, k[line_index], m[line_index])
		p.append(new_point)
	else:
		while True:
			# new space2 when looking at the next line...
			# space2 == space2 - left_r2
			space2 = space2 - r2_left
			# move index to next line
			line_index = line_index + 1
			# move to next line distance. Note that r2_left is the whole
			# line length since previous point was on the previous line
			r2_left = r2[line_index]
			# checks if the point should be in the next line
			if space2 < r2_left:
				# ... then the next point should be in this line
				# finds the point where the distance is space2 
				new_point = calc_next(cs[2*line_index, :], space2, k[line_index], m[line_index])
				p.append(new_point)	
				break		
				
# plotting
plt.figure(figsize=(10,5))
plt.plot(cs[:, 0], cs[:, 1], '-o', color='C0')
pp = np.array(p)
plt.plot(pp[:, 0], pp[:, 1], 'o', color='C1')
plt.grid()
plt.axis('equal')
plt.show()

