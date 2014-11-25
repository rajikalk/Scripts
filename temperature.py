#creates temperature slice plots

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def _Temperature(field, data):
    #This find the temperature, assuming the gas is idea (ie gamma = 5/3)
    #Thermal energy is given in erg/g,
    #Density is given in g/cm^3
    #Number density is given in 1/cm^3
    gamma = data.pf['Gamma']
    Boltzmann_constant = 1.3806488*(10**(-16)) #erg/K
    top= data["ThermalEnergy"] * data["Density"] * (gamma - 1.)
    bottom = data["NumberDensity"] * Boltzmann_constant
    temperature = top/bottom
    inf = np.all(np.isfinite(data["NumberDensity"]))
    return temperature
    
add_field("Temperature", function=_Temperature, units=r"K")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD00*/CE00*.hierarchy")

#save directory
#save_directory = '/media/DATA/YT_output/run1.e-6lessetot_Gcorr_0.75k/xz-plane-Temperature/'
#f = open('temperature_timeseries.txt', 'r+')
#f.write('Time, Average_Temperature \n')
it = 0

#create plots
for pf in ts:
    p = SlicePlot(pf, "z", "Temperature")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    #p.annotate_contour("Density", ncont=9, clim=(1.e-8, 1.e-4), label=True)
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    p.set_zlim('all', zmin = 1e2, zmax = 1e7)
    filename = "frame_" + str(it) + "_time_" + time + ".png"
    p.save(filename)
    #sp = pf.h.sphere(pf.domain_center, (70, 'rsun'))
    it = it + 1
