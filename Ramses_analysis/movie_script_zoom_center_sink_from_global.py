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
import os

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
sink_id = sys.argv[3]
    
files = sorted(glob.glob(input_dir+"data/*/info*.txt"))

#Get zoom_center from input
snkfiles = sorted(glob.glob(input_dir+"data/*/*.snktxt"))
#Find where target sink appears
thres_len = int(sink_id) + 3
for snkfile in snkfiles:
    num_lines = sum(1 for _ in open(snkfile))
    if num_lines > thres_len:
        break
start_file = glob.glob(snkfile.split("stars_output.snktxt")[0]+'info*.txt')[0]
files = files[files.index(start_file):]

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
    
#make projections
for fn in yt.parallel_objects(files):
    fit = files.index(fn)
    save_name = save_dir + "file_frame_" + ("%06d" % fit)+".jpg"
    if os.path.isfile(save_name) == False:
        try:
            ds = yt.load(fn, units_override=units_override)
            
            
            #Get center particle
            sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
            sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
            sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
            sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
            
            
            
            import pdb
            pdb.set_trace()
            region = ds.box(left_corner, right_corner)
            
            time_val = ds.current_time.in_units('yr')
            
            Part_in_region = np.where((region['sink_particle_posx'].in_units('au')>particle_search_bounds_left[0])&(region['sink_particle_posx'].in_units('au')<particle_search_bounds_right[0])&(region['sink_particle_posy'].in_units('au')>particle_search_bounds_left[1])&(region['sink_particle_posy'].in_units('au')<particle_search_bounds_right[1])&(region['sink_particle_posz'].in_units('au')>particle_search_bounds_left[2])&(region['sink_particle_posz'].in_units('au')<particle_search_bounds_right[2]))[0]
            if len(Part_in_region)>0:
                print("DEBUG FROM FILE", fn)
                import pdb
                pdb.set_trace()
            
            proj = yt.ProjectionPlot(ds, 2, ('gas', 'Density'), data_source=region, method='integrate')
            proj.set_buff_size([resolution, resolution])
            proj_array = np.array((proj.frb.data[('gas', 'Density')]/x_width.in_units('cm')).in_units('g/cm**3'))
            
            plt.clf()
            fig, ax = plt.subplots()
            ax.set_xlabel("x (AU)", labelpad=-1, fontsize=12)
            ax.set_ylabel("y (AU)", fontsize=12) #, labelpad=-20
            ax.set_xlim([np.min(X), np.max(X)])
            ax.set_ylim([np.min(Y), np.max(Y)])
            
            cmap=plt.cm.gist_heat
            plot = ax.pcolormesh(X, Y, proj_array, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar(plot, pad=0.0)
            if len(Part_in_region)>0:
                ax.scatter(sink_particle_posx, sink_particle_posy, color='c', s=0.5)
                ax.scatter(sink_particle_posx[most_recent_sink_pos], sink_particle_posy[most_recent_sink_pos], color='g', s=1)
            save_name = save_dir + "file_frame_" + ("%06d" % fit)+".jpg"
            plt.savefig(save_name, format='jpg', bbox_inches='tight')
            print('Created frame on rank', rank, 'at time of', str(time_val), 'to save_dir:', save_name + '.jpg')
        except:
            print("couldn't open file", fn)
    else:
        print(save_name, "already exists")
