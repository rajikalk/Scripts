#Make entropy slices

from yt.mods import *
import numpy as np

def _Entropy(field, data):
    #everything is calculated in cgs.
    Msun = 1.9891*(10**33.) #grams
    m_h = 1.6733e-24        #grams
    k = 1.3806488e-16       #erg/K
    h = 6.6260755e-27       #erg/s
    m = data["CellMassMsun"]*Msun #grams
    N = m/(1.1*m_h) #a number
    V = data["CellVolume"]  #cm^3
    U = data["ThermalEnergy"] * m   #ergs
    S = (N*k)*(np.log((V/N)*(((4*np.pi*m*U)/(3*N*(h**2)))**(3./2.))) + 5./2.)
    inf = np.all(np.isfinite(data["NumberDensity"]))
    return S
    
add_field("Entropy", function=_Entropy, units=r"erg/K")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/hot3/DD00*/CE00*.hierarchy")

#save directory
save_directory = '/media/DATA/YT_output/hot3/xz-plane-Entropy/'
it = 0

#create plots
for pf in ts:
    p = SlicePlot(pf, "y", "Entropy")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    p.set_zlim('all', zmin = 1e32, zmax = 1e38)
    filename = save_directory + "frame_" + str(it) + "_time_" + time + ".png"
    p.save(filename)
    it = it + 1

