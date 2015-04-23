# Simply creates a mass profile

from yt.mods import *
import matplotlib.pyplot as plt

#load data:
pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
time = float(pf.current_time)
#dd = pf.h.all_data()

#Set radius
Msun = 1.9891*(10.**33.) #grams/Msun
Rsun = 6.955e10 #cm in rsun
solar_rad = 140
rad = Rsun*solar_rad #solar radius in centimeters
#sp = pf.sphere("c", (70.0, 'rsun'))

#Get radial velocity
#plot = ProfilePlot(sp, "radius", "density", weight_field="cell_mass")
#plot.set_unit('radius', 'rsun')
# Create a sphere of radius 100 kpc in the center of the box.
my_sphere = pf.h.sphere("c", (50.0, "rsun"))

# Create a profile of the average density vs. radius.
plot = ProfilePlot(my_sphere, "Radius", "Density",
                      weight_field="Density")

# Change the units of the radius into kpc (and not the default in cgs)
plot.set_unit('Radius', 'rsun')
plot.set_log('Radius', False)

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
#save_directory = '/home/reggie/DATA/YT_output/run1.e-6lessetot_Gcorr_0.75k/Density_Profiles/'
filename = 'density_profile_' + str(time) + '.png'
plot.save(filename)
