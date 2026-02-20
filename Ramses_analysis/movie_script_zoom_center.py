#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import my_ramses_module as mym
from mpi4py.MPI import COMM_WORLD as CW
import my_ramses_fields_short as myf
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
    
files = sorted(glob.glob(input_dir+"data/*/info*.txt"))

#Get zoom_center from input
center_pos = [np.nan, np.nan, np.nan]
input_file = input_dir + 'input.nml'
with open(input_file, 'r') as f:
    reader = csv.reader(f, delimiter=" ")
    for row in reader:
        if len(row)>0:
            if row[0] == 'x_center':
                center_pos[0] = float(row[-1])
            if row[0] == 'y_center':
                center_pos[1] = float(row[-1])
            if row[0] == 'z_center':
                center_pos[2] = float(row[-1])
center_pos = yt.YTArray(center_pos, 'pc').in_units('au')
x_width = yt.YTQuantity(20000, 'au')
x = np.linspace(-1*x_width/2, x_width/2, 800)
X, Y = np.meshgrid(x, x)
left_corner = yt.YTArray([center_pos[0].value-(0.5*x_width.value), center_pos[1].value-(0.5*x_width.value), center_pos[2].value-(0.5*x_width.value)], 'AU')
particle_search_bounds_left= []
for left_ind in left_corner:
    if left_ind > 0:
        particle_search_bounds_left.append(left_ind)
    else:
        append_val = yt.YTQuantity(4, 'pc').in_units('au')+left_ind
        particle_search_bounds_left.append(append_val)
right_corner = yt.YTArray([center_pos[0].value+(0.5*x_width.value), center_pos[1].value+(0.5*x_width.value), center_pos[2].value+(0.5*x_width.value)], 'AU')
particle_search_bounds_right= []
for right_ind in right_corner:
    if right_ind > 0:
        particle_search_bounds_right.append(right_ind)
    else:
        append_val = yt.YTQuantity(4, 'pc').in_units('au')+right_ind
        particle_search_bounds_right.append(append_val)

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
    
#make projections
for fn in yt.parallel_objects(files):
    fit = files.index(fn)
    ds = yt.load(fn, units_override=units_override)
    region = ds.box(left_corner, right_corner)
    
    time_val = ds.current_time.in_units('yr')
    
    Part_in_region = np.where((region['sink_particle_posx'].in_units('au')>particle_search_bounds_left[0])&(region['sink_particle_posx'].in_units('au')<particle_search_bounds_right[0])&(region['sink_particle_posy'].in_units('au')>particle_search_bounds_left[1])&(region['sink_particle_posy'].in_units('au')<particle_search_bounds_right[1])&(region['sink_particle_posz'].in_units('au')>particle_search_bounds_left[2])&(region['sink_particle_posz'].in_units('au')<particle_search_bounds_right[2]))[0]
    if len(Part_in_region)>0:
        print("DEBUG FROM FILE", fn)
        import pdb
        pdb.set_trace()
    
    proj = yt.ProjectionPlot(ds, 2, ('gas', 'Density'), data_source=region, method='integrate')
    proj.set_buff_size([800, 800])
    proj_array = np.array((proj.frb.data[('gas', 'Density')]/x_width.in_units('cm')).in_units('g/cm**3'))
    
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xlabel("x (AU)", labelpad=-1, fontsize=12)
    ax.set_ylabel("y (AU)", fontsize=12) #, labelpad=-20
    ax.set_xlim([-1*x_width/2, x_width/2])
    ax.set_ylim([-1*x_width/2, x_width/2])
    
    cmap=plt.cm.gist_heat
    plot = ax.pcolormesh(X, Y, proj_array, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
    plt.gca().set_aspect('equal')
    cbar = plt.colorbar(plot, pad=0.0)
    if len(Part_in_region)>0:
        ax.scatter(sink_particle_posx, sink_particle_posy, color='c', s=0.5)
        ax.scatter(sink_particle_posx[most_recent_sink_pos], sink_particle_posy[most_recent_sink_pos], color='g', s=1)
    save_name = save_dir + "file_frame_" + ("%06d" % fit)
    plt.savefig(save_name + ".jpg", format='jpg', bbox_inches='tight')
    print('Created frame on rank', rank, 'at time of', str(time_val), 'to save_dir:', save_name + '.jpg')
