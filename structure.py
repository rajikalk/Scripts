# mass in orbital plane

from yt.mods import *
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np

#Import all the timesteps for the series:
pf = load("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD0000/CE0000")
ts = TimeSeriesData.from_filenames("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD00*/CE00*.hierarchy")

#Define constants:
dim = pf.domain_dimensions[0]
Msun = pf['MassUnits']       #grams/Msun
lu = pf['LengthUnits']     #length units in cm
tu = pf['TimeUnits']           #time units in sec
DensityUnits = pf['DensityUnits']
gl = lu/dim             #cm in grid length

#Define arrays
time = []   #time
peak = 0.0
TM_1 = []   #total mass in the regions
TM_2 = []
TM_3 = []
TM = []     #total mass in entire grid

#find region bounds
center = 0.3 #Thickness of central slab in grid units
lower_upper_bound = 0.5 - center/2.
upper_lower_bound = 0.5 + center/2.
lower_center = lower_upper_bound/2.
upper_center = 1. - (upper_lower_bound/2.)

it = 0
for pf in ts:
	if it < 100:
	    dd = pf.h.all_data()
	    #Get current time:
	    time.append(pf.current_time*pf["years"])
	    #Define regions in orbital plane and above the plane, and entire grid:
	    grid = pf.h.region([0.5,0.5,0.5],[0.0,0.0,0.0],[1.0,1.0,1.0])
	    reg1 = pf.h.region([0.5,0.5,upper_center],[0.0,0.0,upper_lower_bound],[1.0,1.0,1.0])
	    reg2 = pf.h.region([0.5,0.5,0.5],[0.0,0.0,lower_upper_bound],[1.0,1.0,upper_lower_bound])
	    reg3 = pf.h.region([0.5,0.5,lower_center],[0.0,0.0,0.0],[1.0,1.0,lower_upper_bound])
	    #Find the total mass in these regions.
	    tm = sum(grid["CellMassMsun"])
	    tm_1 = sum(reg1["CellMassMsun"])
	    tm_2 = sum(reg2["CellMassMsun"])
	    tm_3 = sum(reg3["CellMassMsun"])
	    if tm_2 > peak:
		peak = tm_2
	    #Append to mass arrays
	    TM.append(TM)
	    TM_1.append(tm_1)
	    TM_2.append(tm_2)
	    TM_3.append(tm_3)
	    it = it + 1
#time = np.array(time)
#TM = np.array(TM)
#TM_1 = np.array(TM_1)
#TM_2 = np.array(TM_2)
#TM_3 = np.array(TM_3)

print 'Fallback disc mass = ', peak, 'Msun'

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
f.set_size_inches(6,8)
ax1.plot(time, TM_1)

ax1.axvline(x=0.01, color='k')
ax1.axvline(x=0.04, color='k')
ax1.axvline(x=0.10, color='k')
ax1.axvline(x=0.06, color='b', ls='--')
ax1.axvline(x=0.02, color='r', ls=':')
ax1.axvline(x=0.03, color='r', ls=':')
ax1.axvline(x=0.05, color='r', ls=':')
ax1.axvline(x=0.07, color='r', ls=':')
ax1.axvline(x=0.08, color='r', ls=':')
ax1.axvline(x=0.09, color='r', ls=':')

ax1.set_ylabel('Mass ($M_\odot$)')
ax1.set_xlim([0.0, 0.12])
ax1.set_ylim([0.0, 0.3000000001])
ax2.plot(time, TM_2)

ax2.axvline(x=0.01, color='k')
ax2.axvline(x=0.04, color='k')
ax2.axvline(x=0.10, color='k')
ax2.axvline(x=0.06, color='b', ls='--')
ax2.axvline(x=0.02, color='r', ls=':')
ax2.axvline(x=0.03, color='r', ls=':')
ax2.axvline(x=0.05, color='r', ls=':')
ax2.axvline(x=0.07, color='r', ls=':')
ax2.axvline(x=0.08, color='r', ls=':')
ax2.axvline(x=0.09, color='r', ls=':')

ax2.set_ylabel('Mass ($M_\odot$)')
ax3.plot(time, TM_3)

ax3.axvline(x=0.01, color='k')
ax3.axvline(x=0.04, color='k')
ax3.axvline(x=0.10, color='k')
ax3.axvline(x=0.06, color='b', ls='--')
ax3.axvline(x=0.02, color='r', ls=':')
ax3.axvline(x=0.03, color='r', ls=':')
ax3.axvline(x=0.05, color='r', ls=':')
ax3.axvline(x=0.07, color='r', ls=':')
ax3.axvline(x=0.08, color='r', ls=':')
ax3.axvline(x=0.09, color='r', ls=':')

ax3.set_ylabel('Mass ($M_\odot$)')
ax3.set_xlabel('Time ($years$)')
f.subplots_adjust(hspace=0.1)
'''
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.

plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

print 'making plot.'

plt.clf()
fig = plt.figure(1)
print 'made figure.'
fig.set_size_inches(7,9)
print 'set size.'
#fig.yscale("log")
full = fig.add_subplot(411)
print 'made first subplot'
full.plot(time, TM)
print 'plotted first subplot'
#full.set_xlim([0.0, 0.1])
full.set_ylabel('Total Mass ($M_\odot$)')
print 'labelled first subplot'
top = fig.add_subplot(311)
print 'made second subplot'
top.plot(time, TM_1)
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
print 'plotted second subplot'
#top.set_ylim([0.0, 0.15])
#top.set_title('Above orbital plane ('+str(upper_lower_bound)+'>z>1.0)')
#top.set_ylabel('Total Mass ('+str(upper_lower_bound)+'>z>1.0) ($M_\odot$)')
print 'labelled second subplot'
middle = fig.add_subplot(312, sharex=top)
print 'made third subplot'
middle.plot(time, TM_2)
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
print 'plotted third subplot'
#middle.set_title('About orbital plane ('+str(lower_upper_bound)+'>z>'+str(upper_lower_bound)+')')
middle.set_ylabel('Total Mass ($m_\odot$)')
print 'labelled third subplot'
bot = fig.add_subplot(313, sharex=top)
print 'made fourth subplot'
bot.plot(time, TM_3)
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
print 'made fourth subplot'
#bot.set_title('Below orbital plane (0.0>z>'+str(lower_upper_bound)+')')
bot.set_xlabel('time ($years$)')
#bot.set_ylabel('Total Mass (0.0>z>'+str(lower_upper_bound)+') ($M_\odot$)')
print 'labelled fourth subplot'
#fig.tight_layout()
print 'plot figure'
fig.subplots_adjust(hspace = 0)
'''
filename = 'split_mass_profile_30.png'
plt.savefig(filename)
