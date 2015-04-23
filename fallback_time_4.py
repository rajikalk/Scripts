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
particles = pf.h.find_particles_by_type(11)
core_mass = (particles["particle_mass"][0]+particles["particle_mass"][1])

#set array
time = []
distance = []

#Set radius
solar_rad = 10
Rsun = 6.955e10 #cm in Rsun
radius = Rsun*solar_rad #solar radius in centimeters

#calculate gravitational acceleration:
sphere = pf.h.sphere(pf.domain_center, (radius, "cm"))

#Get radial velocity
bin_num = 200
rp = BinnedProfile1D(dd, bin_num, "Radius", 0.0, radius, log_space = False)
rp.add_fields("RadialVelocityABS") #YT gives this in cm/s

i = 0
while i < (bin_num+1):
    #Mass enclosed in this sphere:
    if rp["Radius"][i] < 1.e11:
        rad = 1.e11
    else:
        rad = rp["Radius"][i]
    print "Radius:", rad, "cm"
    sp = pf.h.sphere(pf.domain_center, (rad, "cm"))
    mass = sp.quantities["TotalMass"]()*Msun
    #mass = core_mass
    
    #Calculate gravitational constant
    init_g = g_accel(mass, rad)

    #calculate fallback time:
    rv = rp["RadialVelocityABS"][i]
    init_t = ((-2. * rv)/init_g) #To give time in days
    
    #calculate distance of fallback distance with initial g:
    init_dist = (rv * (init_t/2.)) + (0.5*init_g*((init_t/2.)**2)) + rad
    
    #calculate minimum gravitational acceleration at the peak distance:
    min_g = g_accel(mass, init_dist)
    
    #Find average g:
    ave_g = (init_g+min_g)/2.
    
    #recalculate time:
    t = ((-2. * rv)/ave_g)/(86400.*365.25) #To give time in days
    time.append(t)
    
    #calculate fallback distance with averaged g:
    dist = (((rv * (t/2.)) + (0.5*ave_g*((t/2.)**2.)))+rad)/Rsun
    distance.append(dist)
    
    i = i + 1

#Plot figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(5,5)
#fig.tight_layout()
top = fig.add_subplot(211)
top.plot(rp["Radius"]/length_unit, time)
top.set_title('Fallback time')
top.set_ylabel('Time ($years$)')
bot = fig.add_subplot(212, sharex=top)
bot.plot(rp["Radius"]/length_unit, distance)
bot.set_title('Distance')
bot.set_ylabel('Distance ($R_\odot$)')
bot.set_xlabel('Radius ($R_\odot$)')
filename = 'rotating_time '+str(pf.current_time)+' fallback_time_profile.png'
plt.savefig(filename)

