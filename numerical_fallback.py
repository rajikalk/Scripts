#Numeral simulation of fallback
#Last editted 1/5/2014

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
#Find the change in height for a projectile over a time step. All units given in CGS
    first_term = velocity*time_step
    second_term = 0.5*acceleration*(time_step**2.)
    h = first_term + second_term
    return h
    
def velocity(initial_velocity, acceleration, time_step):
#Calculates the velocity over a given time step. All units given in CGS.
    v = initial_velocity + acceleration*time_step
    return v

#load data:
pf = load("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD0013/CE0013")
dd = pf.h.all_data()

#Get the grids in the data dump:
grids = pf.h.grids[0]

#Define Arrays:
radius = []
time = []
peak_r = []
    
#system:
Msun = 1.9891*(10.**33.) #grams/Msun
Rsun = 6.955e10 #cm in Rsun
pmass = (dd['ParticleMassMsun'][0] + dd['ParticleMassMsun'][1])*Msun
dt = 60. #Time step
day = 60.*60.*24. #seconds in a day

#Get radial profile velocity
half_width = (pf.domain_width[0]*pf['rsun'])/2.0
diagonal = (half_width**2.0 + half_width**2.0 + half_width**2.0)**0.5
bin_num = 200
rad = diagonal*Rsun
rp = BinnedProfile1D(dd, bin_num, "Radius", 0.0, rad, log_space = False)
rp.add_fields("RadialVelocityABS") #YT gives this in cm/s
rp.add_fields("RadialVelocityKMS")
i=0

while i < (bin_num+1): #Iterates through all bins.
    print "For i = ", i
    init_r = rp["Radius"][i] #Takes initial radius from binned profile.
    radius.append(init_r/Rsun) #Adds this too the radius array
    r = init_r #Sets radius
    v_rad = rp["RadialVelocityABS"][i] #Takes initial radial velocity from binned profile.
    a = g_accel(pmass, r)  #Works out initial acceleration due to gravity. Assumes gas is negligible.
    t = 0. #time passed
    final_dr = (-1.)*r #Final radius should be at the centre, hence why it is negative initial r.
    dr = 0.
    peak = 0.
    while r > 0.0:
        height = projectile_height(v_rad, dt, a)
        r = r + height
        if r > peak:
            peak = r
        dr = r - init_r
        v_rad = velocity(v_rad, a, dt)
        a = g_accel(pmass, r)
        t = t + dt
    time.append(t/day)
    peak_r.append(peak/Rsun)
    i = i + 1
print "radius size:", len(radius), "time size:", len(time), "peak size:", len(peak_r)

fig = plt.figure(1)
fig.clf()
#fig.set_size_inches(7,9)
#fig.yscale("log")
top = fig.add_subplot(211)
top.plot(radius, time)
#top.set_xlim([0.0, 100.0])
top.set_title('Fallback time of gas at that radius')
top.set_ylabel('Time (days)')
bot = fig.add_subplot(212, sharex=top)
bot.plot(radius, peak_r)
bot.set_title('Peak height gas with reach')
bot.set_xlabel('Radius ($R_\odot$)')
bot.set_ylabel('Peak height ($R_\odot$)')
fig.tight_layout()
filename = '/media/DATA/YT_output/run1.e-6lessetot_Gcorr_0.75k/fallback_time_profile.png'
plt.savefig(filename)


