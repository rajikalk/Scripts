#Fallback time along one line of sight
from yt.mods import *
import matplotlib.pyplot as plt

#Defines function:
def g_accel(mass, radius):
#Takes mass and radius in in cgs adn spits out value in cgs
    g_top = (-6.67259e-8)*mass
    g_bot = (radius)**2.
    g = g_top/g_bot
    return g

#load data:
pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr/DD0013/CE0013")
dd = pf.h.all_data()

#Get the grids in the data dump:
g = pf.h.grids[0]
    
#system:
Msun = 1.9891*(10**33) #grams/Msun
Rsun = 6.955e10 #cm in Rsun
pmass = (dd['ParticleMassMsun'][0] + dd['ParticleMassMsun'][1])*Msun

#set inital radius in grid units:
z = 0.51

#Define arrays:
distance = []
change = []
time =[]
radius = []

while z < 1.0:
    #set radius:
    radius_cm = (z - 0.5)*pf['cm']
    
    #Mass enclosed in this sphere:
    if radius_cm < 1.e11:
        rad = 1.e11
    else:
        rad = radius_cm
    print "Radius:", rad, "cm"
    radius_Rsun = rad/Rsun
    radius.append(radius_Rsun)
    sp = pf.h.sphere(pf.domain_center, (rad, "cm"))
    mass = sp.quantities["TotalMass"]()*Msun + pmass
    
    #Calculate gravitational acceleration:
    init_g = g_accel(mass, rad)

    #calculate fallback time:
    rv = abs(g['z-velocity'][z, 0.5, 0.5])
    print "Radial velocity:", rv, "cm/s"
    init_t = ((-2. * rv)/init_g) #To give time in days
    
    #calculate distance of fallback distance with initial g:
    init_dist = (rv * (init_t/2.)) + (0.5*init_g*((init_t/2.)**2)) + rad
    
    #calculate minimum gravitational acceleration at the peak distance:
    min_g = g_accel(mass, init_dist)
    
    #Find average g:
    ave_g = (init_g+min_g)/2.
    
    #recalculate time:
    t = ((-2. * rv)/ave_g)/(86400.) #To give time in days
    time.append(t)
    
    #calculate fallback distance with averaged g:
    delta_dist = ((rv * (t/2.)) + (0.5*ave_g*((t/2.)**2.)))/Rsun
    change.append(delta_dist)
    #print "Change in height:", delta_dist, "cm"
    dist = rad/Rsun + delta_dist
    distance.append(dist)
    
    z = z + 0.01
    
#Plot figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(9,12)
#fig.tight_layout()
top = fig.add_subplot(311)
top.plot(radius, time)
top.set_title('Fallback time')
top.set_ylabel('Time ($days$)')
mid = fig.add_subplot(312, sharex=top)
mid.plot(radius, distance)
mid.set_ylabel('Distance ($R_\odot$)')
bot = fig.add_subplot(313, sharex=top)
bot.plot(radius, change)
bot.set_ylabel('$\Delta$Distance ($R_\odot$)')
bot.set_xlabel('Radius ($R_\odot$)')
filename = 'method1-2-z_nonrotating_time '+str(pf.current_time)+' fallback_time_profile.png'
plt.savefig(filename)


