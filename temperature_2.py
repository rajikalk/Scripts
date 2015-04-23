#creates temperature slice plots

from yt.mods import *
import numpy as np

def _Temperature(field, data):
    #This find the temperature, assuming the gas is idea (ie gamma = 5/3)
    #Thermal energy is given in erg/g,
    #Density is given in g/cm^3
    #Number density is given in 1/cm^3
    Msun = 1.9891*(10**33.) #grams in the sun
    gamma = 5./3.
    Boltzmann_constant = 1.3806488*(10**(-16)) #erg/K
    top= data["ThermalEnergy"] * data["Density"] * (gamma - 1.)
    bottom = data["NumberDensity"] * Boltzmann_constant
    temperature = top/bottom
    inf = np.all(np.isfinite(data["NumberDensity"]))
    return temperature
    
add_field("Temperature", function=_Temperature, units=r"K")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/hot2/DD00*/CE00*.hierarchy")

#save directory
save_directory = '/media/DATA/YT_output/hot2/xz-plane-Temperature/'
it = 0

#create plots
for pf in ts:
    p = SlicePlot(pf, "y", "Temperature")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    p.set_zlim('all', zmin = 1e2, zmax = 1e7)
    filename = save_directory + "frame_" + str(it) + "_time_" + time + ".png"
    p.save(filename)
    it = it + 1
