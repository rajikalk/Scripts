#Numeral simulation of fallback

from yt.mods import *
import matplotlib.pyplot as plt

#Defines function:
def g_accel(mass, radius):
#Takes mass and radius in cgs and spits out value in cgs
    g_top = (-6.67259e-8)*mass
    g_bot = (radius)**2.
    g = g_top/g_bot
    return g
    
def projectile_height(velocity, time_step, acceleration):
    first_term = velocity*time_step
    second_term = 0.5*acceleration*(time_step**2.)
    h = first_term + second_term
    return h
    
def velocity(initial_velocity, acceleration, time_step):
    v = initial_velocity + acceleration*time_step
    return v

#load data:
pf = load("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD0013/CE0013")
dd = pf.h.all_data()

#Get the grids in the data dump:
grids = pf.h.grids[0]

#Define Arrays:
time = []
distance = []
distance_cm = []
change = []
    
#system:
Msun = 1.9891*(10.**33.) #grams/Msun
Rsun = 6.955e10 #cm in Rsun
pmass = (dd['ParticleMassMsun'][0] + dd['ParticleMassMsun'][1])*Msun
days = 60.*60.*24. #seconds in a day

#Set inital values
init_radius = 0.83*pf['cm']
init_rsun = 0.83*pf['rsun']
radius = init_radius
g = g_accel(pmass, radius)
rv = abs(grids['y-velocity'][0.5, 0.83, 0.5])
print "inital RV:", rv, "cm/s"
dt = 60. #in seconds
change_val = 0.0
i = 0

#Append inital values to arrays:
time.append(0)
distance.append(init_rsun)
distance_cm.append(init_radius)
change.append(change_val)

while radius >= 0.0:
    time_val = (time[i] + dt)
    time.append(time_val)
    dh = projectile_height(rv, dt, g)
    radius = distance_cm[i] + dh
    distance_cm.append(radius)
    radius_rsun = radius/Rsun
    if radius_rsun >= 0.0:
        distance.append(radius_rsun)
    else:
        distance.append(0)
    change_val = radius - init_radius
    change_val_rsun = change_val/Rsun
    change.append(change_val_rsun)
    
    g = g_accel(pmass, radius)
    rv = velocity(rv, g, dt)
    
    print "time:", time_val, "sec"
    print "radius", radius, "cm"
    print "change", change_val, "cm"
    
    i = i + 1
    
for i in range(len(time)):
    time[i] = time[i]/days

#Plot figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(9,9)
#fig.tight_layout()
top = fig.add_subplot(211)
top.plot(time, distance)
top.set_title('Distance of fallback')
top.set_ylabel('Distance ($R_\odot$)')
bot = fig.add_subplot(212, sharex=top)
bot.plot(time, change)
bot.set_ylabel('$\Delta$Distance ($R_\odot$)')
bot.set_xlabel('Time ($days$)')
bot.set_ylim(0.0, 2.0)
filename = '/media/DATA/YT_output/run1.e-6lessetot_Gcorr_0.75k/fallback_profile.png'
plt.savefig(filename)

