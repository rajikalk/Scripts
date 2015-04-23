from yt.mods import *
import matplotlib.pyplot as plt

#Load DataCube
pf = load("/media/DATA/Simulations/smallbox/rotation/run1.e-6/DD0025/CE0025")

#system:
length_unit = 6.955e10 #1Rsun in cm
system_radius = 6.955e10*100 #4au in cm
#int_sphere = pf.h.sphere(pf.domain_center, (4, "au"))
#bulk_vel = int_sphere.quantities["BulkVelocity"]()
sphere = pf.h.sphere(pf.domain_center, (4, "au"))
#sphere.set_field_parameter("bulk_velocity", bulk_vel)

#create radial profile:
radial_profile = BinnedProfile1D(sphere, 200, "Radius", 0.0, system_radius, log_space=False)
radial_profile.add_fields("RadialVelocityABS")

plt.figure(1) #working on display 1
plt.clf()
#plt.ylim(ymin=-1.0e-1,ymax=1.0e7) #sets the limits of the y axis
#plt.xlim(xmin=1.0e1,xmax=2.0e2) #sets the limits of the x axis
plt.xscale("log")
plt.yscale("log")
plt.title("Radial Velocity vs Radius (t="+str(pf.current_time)+" $yr$)")
plt.xlabel("Radius ($R_{\odot}$)")
plt.ylabel("Radial Velocity ($ms^{-1}$)")
plt.plot(radial_profile["Radius"]/length_unit,radial_profile["RadialVelocityABS"],'-b') 
filename = 'time '+str(pf.current_time)+' radial_velocity_profile.png'
plt.savefig(filename)


