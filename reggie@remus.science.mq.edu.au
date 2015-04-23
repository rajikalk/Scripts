#Creates projection plots of each data dump.

from yt.mods import *
import matplotlib.pyplot as plt

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD00*/CE00*.hierarchy")

#save_directory = '~/YT_output/run1.e-6lessetot_Gcorr_0.75k/xy-plane/'
it = 0

for pf in ts:
    p = SlicePlot(pf, "z", "Density")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    time = str(float(pf.current_time))
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    p.set_zlim('all', zmin = 1e-9, zmax = 1e-4)
    filename = "frame_" + ("%03d" % it) + ".png"
    p.save(filename)
    it = it + 1
