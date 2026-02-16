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
m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=5000., start_time=-500)

usable_files = mym.find_files(m_times, files)
sink_form_time = mym.find_sink_formation_time(files)
del files

if os.path.isfile(save_dir + output) == False:
    f = open(save_dir + output, 'a+')
    f.close()
f = open(save_dir + output, 'w')
f.write('Time, Particles, Gas\n')
f.close()
storage = {}

ts = yt.DatasetSeries(usable_files, parallel=True)
for sto, ds in ts.piter(storage=storage):
    time = ds.current_time.in_units('yr').value - sink_form_time
    dd = ds.all_data()
    center_pos = dd['Center_Position']
    center_vel = dd['Center_Velocity']
    part_pos = dd['All_Particle_Positions']
    part_mass = dd['All_Particle_Masses']
    part_vel = dd['All_Particle_Velocities']
    print("part_mass =", myf.get_part_mass())
    print("part_vel =", myf.get_part_vel())
    print("center_vel =", myf.get_center_vel())
    print("part_pos =", myf.get_part_pos())
    print("center_pos =", myf.get_center_pos())
    if ('io', 'particle_mass') in ds.derived_field_list:
        '''
        F_rad = dd['Gravitational_Force_on_particles_Rad']
        F_tan = dd['Gravitational_Force_on_particles_Tan']
        for part in range(len(F_rad)):
            write_data.append(F_rad[part].value)
            write_data.append(F_tan[part].value)
        '''
        part_ang = np.sum(dd['Particle_Specific_Angular_Momentum'].value)
    else:
        part_ang = 0.0
    gas_ang = np.sum(dd['Specific_Angular_Momentum'].value)
    write_data = [time, part_ang, gas_ang]

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
    for key, value in storage.items():
        if value is None:
            del_keys.append(key)
    for key in del_keys:
        del storage[key]
    f = open(save_dir + output, 'a')
    for it in sorted(np.array(list(storage.keys())).astype(np.float)):
        f.write(storage[str(it)])
        print("Printed line:", storage[str(it)])
    f.close()
print("finished saving data to", save_dir + output)
