#Fallback time
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
pf = load("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD0013/CE0013")
dd = pf.h.all_data()

#system:
length_unit = 6.955e10 #cm in rsun
Msun = 1.9891*(10**33) #grams/Msun
mass = (dd['ParticleMassMsun'][0] + dd['ParticleMassMsun'][1])*Msun
days = 60.*60.*24. #seconds in a day

#set array
time = []
distance = []
change = []

#Set radius
solar_rad = 140
Rsun = 6.955e10 #cm in Rsun
radius = Rsun*solar_rad #solar radius in centimeters

#calculate gravitational acceleration:
sphere = pf.h.sphere(pf.domain_center, (radius, "cm"))

#Get radial velocity
bin_num = 200
rp = BinnedProfile1D(dd, bin_num, "Radius", 0.0, radius, log_space = False)
rp.add_fields("RadialVelocity") #YT gives this in cm/s

i = 0
while i < (bin_num+1):
    #Mass enclosed in this sphere:
    if rp["Radius"][i] < 1.e11:
        rad = 1.e11
    else:
        rad = rp["Radius"][i]
    print "Radius:", rad, "cm"
    sp = pf.h.sphere(pf.domain_center, (rad, "cm"))
    
    #Calculate gravitational acceleration:
    g = g_accel(mass, rad)

    #calculate fallback time:
    rv = rp["RadialVelocity"][i]
    
    #calculate possible fallback times
    t1 = (-rv + ((rv**2.)- 2.*g*rad)**0.5)/g
    t2 = (-rv - ((rv**2.)- 2.*g*rad)**0.5)/g
    
    if t1 >= 0:
        t = t1
    else:
        t = t2
    t_days = t/days
    time.append(t_days)
    
    #calculate furthest distance:
    h = (-rv)/(2.*g)
    change.append(h)
    
    height = h + rad/Rsun
    distance.append(height)
    
    i = i + 1

#Plot figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(9,12)
#fig.tight_layout()
top = fig.add_subplot(311)
top.plot(rp["Radius"]/length_unit, time)
top.set_title('Fallback time')
top.set_ylabel('Time ($days$)')
mid = fig.add_subplot(312, sharex=top)
mid.plot(rp["Radius"]/length_unit, distance)
mid.set_ylabel('Distance ($R_\odot$)')
bot = fig.add_subplot(313, sharex=top)
bot.plot(rp["Radius"]/length_unit, change)
bot.set_ylabel('$\Delta$Distance ($R_\odot$)')
bot.set_xlabel('Radius ($R_\odot$)')
filename = 'method3_rotating_time '+str(pf.current_time)+' fallback_time_profile.png'
plt.savefig(filename)
