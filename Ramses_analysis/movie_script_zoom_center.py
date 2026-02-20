#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import gc
import my_ramses_fields_short as myf
import csv

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
files = sorted(glob.glob(input_dir+"data/*/info*.txt"))
gc.collect()

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
x_width = 5000
left_corner = yt.YTArray([center_pos[0].value-(0.75*x_width), center_pos[1].value-(0.75*x_width), center_pos[2].value-(0.75*x_width)], 'AU')
right_corner = yt.YTArray([center_pos[0].value+(0.75*x_width), center_pos[1].value+(0.75*x_width), center_pos[2].value+(0.75*x_width)], 'AU')

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})
    
#make projections
for fn in yt.parallel_objects(files):
    ds = yt.load(fn, units_override=units_override)
    region = ds.box(left_corner, right_corner)
    
    import pdb
    pdb.set_trace()
