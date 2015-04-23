#Creates projection plots of each data dump.

from yt.mods import *
import matplotlib.pyplot as plt

#Import all the timesteps for the series:
ts = DatasetSeries("/Users/rajikak/Output/CircumbinaryOutFlow_0.25/WIND_hdf5_plt_cnt_*", particle_filename="/Users/rajikak/Output/CircumbinaryOutFlow_0.25/WIND_hdf5_part_*")

#save_directory = '~/YT_output/CircumbinaryOutFlow/Density_xz_plane/'
it = 440
max_file = 453

while it < max_file:
    it_str = ("%04d" % it)
    pf = load("/Users/rajikak/Output/CircumbinaryOutFlow_0.25/WIND_hdf5_plt_cnt_"+it_str, particle_filename="/Users/rajikak/Output/CircumbinaryOutFlow_0.25/WIND_hdf5_part_"+it_str)
    p = ProjectionPlot(pf, "z", "density")
    p.set_log("density", False)
    p.zoom(10)
    p.annotate_streamlines('magx', 'magy')
    p.annotate_velocity(normalize=True)
    if len(pf.field_list) > 12:
        p.annotate_particles(10.0, p_size = 50)
    p.set_cmap(field="density", cmap="hot")
    time = str(float(pf.current_time))
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    #p.set_zlim('all', zmin = 1e-20, zmax = 1e-13)
    filename = "frame_" + ("%06d" % it) + ".png"
    p.save(filename)
    it = it + 1