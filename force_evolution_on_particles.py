#!/usr/bin/env python
import sys
import yt
yt.enable_parallelism()
import my_fields as myf
import my_module as mym
import glob
import os
import numpy as np
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="What will you save your output files as?", default="force.csv")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

rank = CW.Get_rank()

# Read in directories:
path = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

args = parse_inputs()
output = args.output
files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=5000.)

usable_files = mym.find_files(m_times, files)
sink_form_time = mym.find_sink_formation_time(files)
del files

if os.path.isfile(save_dir + output) == False:
    f = open(save_dir + output, 'a+')
    f.close()
f = open(save_dir + output, 'w')
f.write('Time, Particle 1: F_rad, F_tan, Particle 2: F_rad, F_tan\n')
f.close()
storage = {}

ts = yt.DatasetSeries(usable_files, parallel=True)
for sto, ds in ts.piter(storage=storage):
    if ('io', 'particle_mass') in ds.derived_field_list:
        time = ds.current_time.in_units('yr').value - sink_form_time
        dd = ds.all_data()
        F_rad = dd['Gravitational_Force_on_particles_Rad']
        F_tan = dd['Gravitational_Force_on_particles_Tan']
        write_data = [time]
        for part in range(len(F_rad)):
            write_data.append(F_rad[part].value)
            write_data.append(F_tan[part].value)
        
        write_string = ''
        for dit in range(len(write_data)):
            if dit != len(write_data)-1:
                write_string = write_string + str(write_data[dit]) + ','
            else:
                write_string = write_string + str(write_data[dit]) + '\n'
        sto.result = write_string
        sto.result_id = str(time)

if rank == 0:
    del_keys = []
    for key, value in storage.iteritems():
        if value is None:
            del_keys.append(key)
    for key in del_keys:
        del storage[key]
    f = open(save_dir + output, 'a')
    for it in sorted(np.array(storage.keys()).astype(np.float)):
        f.write(storage[str(it)])
        print "Printed line:", storage[str(it)]
    f.close()
