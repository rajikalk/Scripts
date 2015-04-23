#Creates a density profile (Based onRoberto's script):

from yt.mods import *
import numpy as np
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt

def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

#Load DataCube
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD00*/CE00*.hierarchy")

#system:
#Msun = init_pf['MassUnits']       #grams/Msun
#lu = init_pf['LengthUnits']     #length units in cm
#tu = init_pf['TimeUnits']           #time units in sec
#DensityUnits = init_pf['DensityUnits']
lu = 10000000000000 #1Rsun in cm
system_radius = 6.0e13 #4au in cm
radius = 4 #in rsun
Rsun = 6.955e10
it = 0
#arrays:
times = []
gradient = []
t_min = 0.
grad_min = 0.

for pf in ts:
	time = float(pf.current_time)
	times.append(time)
	dd = pf.h.all_data()
	center = [float(pf.domain_center[0]), float(pf.domain_center[1]), float(pf.domain_center[2])]
	pp1 = [float(dd['particle_position_x'][0]), float(dd['particle_position_y'][0]), float(dd['particle_position_z'][0])]
	pp2 = [float(dd['particle_position_x'][1]), float(dd['particle_position_y'][1]), float(dd['particle_position_z'][1])]
	pd1 = distance(pp1, center)*(pf['LengthUnits']/Rsun)
	pd2 = distance(pp2, center)*(pf['LengthUnits']/Rsun)
	sep = distance(pp1, pp2)*pf['LengthUnits']

	#create sphere for profile:
	#sp1 = pf.h.sphere(pp1, (sep/2., 'cm'))
	#sp2 = pf.h.sphere(pp2, (sep/2., 'cm'))
	#prof1 = BinnedProfile1D(sp1, 100, "Radius", 0.0, sep/2., log_space=False)
	#prof1.add_fields(["Density"], weight="Density")
	#prof2 = BinnedProfile1D(sp2, 100, "Radius", 0.0, sep/2., log_space=False)
	#prof2.add_fields(["Density"], weight="Density")
	sp = pf.h.sphere(center, (50.0, 'rsun'))
	prof = BinnedProfile1D(sp, 100, "Radius", 0.0, (50.0/pf['rsun'])*pf['LengthUnits'], log_space=False)
	prof.add_fields(["Density"], weight="Density")
	prof["Radius"] = prof["Radius"]/Rsun
	
	#Calculate gradient
	rise = prof["Density"][40] - prof["Density"][0]
	run = prof["Radius"][40] - prof["Radius"][0]
	grad = rise/run
	if grad < grad_min:
		grad_min = grad
		t_min = time
	gradient.append(grad)
	x = [prof["Radius"][0], prof["Radius"][40], (-prof["Density"][0])/grad]
	y = [prof["Density"][0], prof["Density"][40], 0.0]

	print "time:", time, ", grad:", grad, "g/cm^3/Rsun"

	plt.clf()
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.axvline(x=pd1, ls=':')
	plt.axvline(x=pd2, ls=':')
	#plt.plot(prof1['Radius'], prof1['Density'], 'k')
	#plt.plot(prof2['Radius'], prof2['Density'], 'b')
	plt.plot(prof['Radius'], prof['Density'], 'k')
	plt.plot(x, y, 'r--')	
	plt.xlabel("Radius ($R_\odot$)")
	plt.ylabel("Density ($gcm^{-3}$)")
	#plt.yscale('log')
	plt.ylim([0.0, 6.e-5])
	plt.xlim([0.0, 50.0])
	#plt.title("Current Time:" + str(time) + "years")
	filename = 'linear_frame_' + str(it) + '_time_' + str(time) + '.png'
	plt.savefig(filename)
	it = it + 1
	
plt.clf()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(times, gradient)
plt.ylim([-2.e-6, 0.5e-6])
plt.axvline(x=0.01, color='k')
plt.axvline(x=0.04, color='k')
plt.axvline(x=0.10, color='k')
plt.axvline(x=0.06, color='b', ls='--')
plt.axvline(x=0.02, color='r', ls=':')
plt.axvline(x=0.03, color='r', ls=':')
plt.axvline(x=0.05, color='r', ls=':')
plt.axvline(x=0.07, color='r', ls=':')
plt.axvline(x=0.08, color='r', ls=':')
plt.axvline(x=0.09, color='r', ls=':')
plt.axvline(x=0.06, color='b', ls='--')
plt.xlabel("Time ($years$)")
plt.ylabel("Gradient ($gcm^{-3}R_{\odot}^{-1}$)")
plt.savefig("grad_timeseries.png")
print "t_min:", t_min, "years, grad_min:", grad_min

