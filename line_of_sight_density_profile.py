from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/Hot2/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/Hot2/DD0000/CE0000")

dim = pf.domain_dimensions[0]
it = 0

for pf in ts:
	
	#import all data
	dd = pf.h.all_data()
	g = pf.h.grids[0]

	#get particle positions
	pp1 = [dd['particle_position_x'][0]*dim, dd['particle_position_y'][0]*dim, dd['particle_position_z'][0]*dim]
	pp2 = [dd['particle_position_x'][1]*dim, dd['particle_position_y'][1]*dim, dd['particle_position_z'][1]*dim]

	#find particle distance
	pd1 = ((pp1[0]-dim/2.)**2. + (pp1[1]-dim/2.)**2. + (pp1[2]-dim/2.)**2.)**0.5
	pd2 = ((pp2[0]-dim/2.)**2. + (pp2[1]-dim/2.)**2. + (pp2[2]-dim/2.)**2.)**0.5

	#calculate increment vector
	frac1 = 1./(((pp1[0]-dim/2.)**2. + (pp1[1]-dim/2.)**2.)**0.5)
	frac2 = 1./(((pp2[0]-dim/2.)**2. + (pp2[2]-dim/2.)**2.)**0.5)
	sv1 = [(pp1[0]-dim/2.)*frac1, (pp1[1]-dim/2.)*frac1, (pp1[2]-dim/2.)*frac1]
	sv2 = [(pp2[0]-dim/2.)*frac2, (pp2[1]-dim/2.)*frac2, (pp2[2]-dim/2.)*frac2]

	print "Particle 1:"

	dist1 = 0.0
	pos1 = [dim/2., dim/2., dim/2.]
	distance1 = []
	density1 = []
	while dist1 < 128.0:
		den = g['Density'][pos1[0], pos1[1], pos1[2]]
		density1.append(den)
		distance1.append((dist1/dim)*pf['rsun'])

		pos1 = [pos1[0]+sv1[0], pos1[1]+sv1[1], pos1[2]+sv1[2]]
		dist1 = ((pos1[0]-dim/2.)**2. + (pos1[1]-dim/2.)**2. + (pos1[2]-dim/2.)**2.)**0.5
		print dist1

	print "Particle 2:"
	dist2 = 0.0
	pos2 = [dim/2., dim/2., dim/2.]
	distance2 = []
	density2 = []
	while dist2 < 128.0:
		den = g['Density'][pos2[0], pos2[1], pos2[2]]
		density2.append(den)
		distance2.append((dist2/dim)*pf['rsun'])

		pos2 = [pos2[0]+sv2[0], pos2[1]+sv2[1], pos2[2]+sv2[2]]
		dist2 = ((pos2[0]-dim/2.)**2. + (pos2[1]-dim/2.)**2. + (pos2[2]-dim/2.)**2.)**0.5
		print dist2

	time = pf.current_time

	plt.clf()
	plt.plot(distance1, density1, 'k--')
	plt.plot(distance2, density2, 'b--')
	plt.axvline(x=pd1, color='k')
	plt.axvline(x=pd2, color='b')
	plt.xlabel("Radius ($R_\odot$)")
	plt.ylabel("Density ($gcm^{-3}$)")
	plt.savefig("LOS_density_frame_"+str(it)+"_time_"+str(time)+".png")

	it = it + 1
