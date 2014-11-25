#finds density gradient across particles

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

def rel_position(pos1, pos2):
	#finds the relative position of the position to the CoM in grid coordinates as used in the code.
	#this is just to work out the sign of the linear momentum
	x = pos1[0] - pos2[0]
	y = pos1[1] - pos2[1]
	z = pos1[2] - pos2[2]
	rel = [x,y,z]
	return rel

def gradient(y2, y1, x2, x1):
	rise = y2 - y1
	run = x2- x1
	grad = rise/run
	return grad

def distance(point1, point2):
    #Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot3/DD0000/CE0000")
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot3/DD00*/CE00*.hierarchy")

dim = pf.domain_dimensions[0]
center = [dim/2., dim/2., dim/2.]
time = []
den_grad_p1 = []
den_grad_p2 = []

for pf in ts:
	t_val = pf.current_time
	time.append(t_val)
	dd = pf.h.all_data()
	g = pf.h.grids[0]
	    
	#Find the position of the particles:
	pp1 = [dd['particle_position_x'][0]*dim, dd['particle_position_y'][0]*dim, dd['particle_position_z'][0]*dim]
	pp2 = [dd['particle_position_x'][1]*dim, dd['particle_position_y'][1]*dim, dd['particle_position_z'][1]*dim]
	
	#Find relative position:
	rel1 = rel_position(pp1, center)
	rel2 = rel_position(pp2, center)
	
	#Find inner and outer points:
	I1 = [pp1[0] - 0.5*rel1[0], pp1[1] - 0.5*rel1[1], pp1[2] - 0.5*rel1[2]]
	I2 = [pp2[0] - 0.5*rel2[0], pp2[1] - 0.5*rel2[1], pp2[2] - 0.5*rel2[2]]
	O1 = [pp1[0] + 0.5*rel1[0], pp1[1] + 0.5*rel1[1], pp1[2] + 0.5*rel1[2]]
	O2 = [pp2[0] + 0.5*rel2[0], pp2[1] + 0.5*rel2[1], pp2[2] + 0.5*rel2[2]]
	
	#Find distance between points:
	dist1 = ((distance(I1, O1))/dim)*pf['rsun']
	dist2 = ((distance(I2, O2))/dim)*pf['rsun']
	
	#Find density difference:
	den1 = g['Density'][O1[0], O1[1], O1[2]] - g['Density'][I1[0], I1[1], I1[2]]
	den2 = g['Density'][O2[0], O2[1], O2[2]] - g['Density'][I2[0], I2[1], I2[2]]
	
	#Find gradient:
	grad1 = abs(den1/dist1)
	grad2 = abs(den2/dist2)
	den_grad_p1.append(grad1)
	den_grad_p2.append(grad2)
	
	print "Gradient across P1:", grad1
	print "Gradient across P2:", grad2

plt.clf()
plt.plot(time, den_grad_p1, 'k-')
plt.plot(time, den_grad_p2, 'b--')
plt.xlabel("Time ($years$)")
plt.ylabel("Density Gradient across particle")
#plt.yscale('log')
#plt.ylim([1.e-7, 5.e-5])
#plt.title("Current Time:" + str(time) + "years")
filename = 'hot3_grad.png'
plt.savefig(filename)

