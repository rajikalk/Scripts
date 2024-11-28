import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
import glob
import my_flash_module as mym
import argparse

sink_evol_pickle = sys.argv[1]
primary_sink = sys.argv[2]
secondary_sink = sys.argv[3]
sim_dir = '/scratch/ek9/ccf100/sf_outflow/r1024mM5Ma2A1oh/'

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-box", "--box_length", help="how big do you want your resimulate box?", type=float, default=0.1)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()
box_length = yt.YTQuantity(args.box_length, 'pc')

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

#Get binary CoM when secondary have 0.2Msun

secondary_mass_ind = np.argmin(abs(yt.YTArray(sink_data[secondary_sink]['mass'], 'g').in_units('Msun').value - 0.2))
binary_characteristic_time = yt.YTQuantity(sink_data[secondary_sink]['time'][secondary_mass_ind], 's')
binary_masses = yt.YTArray([sink_data[primary_sink]['mass'][secondary_mass_ind], sink_data[secondary_sink]['mass'][secondary_mass_ind]], 'g')
binary_positions = yt.YTArray([[sink_data[primary_sink]['posx'][secondary_mass_ind], sink_data[secondary_sink]['posx'][secondary_mass_ind]], [sink_data[primary_sink]['posy'][secondary_mass_ind], sink_data[secondary_sink]['posy'][secondary_mass_ind]], [sink_data[primary_sink]['posz'][secondary_mass_ind], sink_data[secondary_sink]['posz'][secondary_mass_ind]]], 'cm')
binary_velocities = yt.YTArray([[sink_data[primary_sink]['velx'][secondary_mass_ind], sink_data[secondary_sink]['velx'][secondary_mass_ind]], [sink_data[primary_sink]['vely'][secondary_mass_ind], sink_data[secondary_sink]['vely'][secondary_mass_ind]], [sink_data[primary_sink]['velz'][secondary_mass_ind], sink_data[secondary_sink]['velz'][secondary_mass_ind]]], 'cm/s')
CoM_pos = (binary_positions.T[0] * binary_masses[0] + binary_positions.T[1] * binary_masses[1])/np.sum(binary_masses)
CoM_vel = (binary_velocities.T[0] * binary_masses[0] + binary_velocities.T[1] * binary_masses[1])/np.sum(binary_masses)


#Find frame before primary formation
form_time = yt.YTArray(sink_data[primary_sink]['time'][0], 's')
import pdb
pdb.set_trace()

'''
xmin = form_position[0] - box_length/2
xmax = form_position[0] + box_length/2
ymin = form_position[1] - box_length/2
ymax = form_position[1] + box_length/2
zmin = form_position[2] - box_length/2
zmax = form_position[2] + box_length/2

#Find simfile right before sink forms
files = sorted(glob.glob(sim_dir + '*plt_cnt*'))
files = [ x for x in files if "_proj_" not in x ]
form_file = mym.find_files([form_time.value], files)
prev_ind = int(form_file[0].split('_')[-1]) - 1
prev_file = form_file[0][:-4] + str(prev_ind)

#Calculate bulk file
part_file = 'part'.join(prev_file.split('plt_cnt'))
ds = yt.load(prev_file, particle_filename=part_file)
left_corner = yt.YTArray([xmin, ymin, zmin], 'cm')
right_corner = yt.YTArray([xmax, ymax, zmax], 'cm')
region = ds.box(left_corner, right_corner)

sim_velx_offset = -1 * np.mean(region['velx'])
sim_vely_offset = -1 * np.mean(region['vely'])
sim_velz_offset = -1 * np.mean(region['velz'])

#Print results
print("sim_input_file = ", prev_file)
print("xmin = " + str(xmin))
print("xmax = " + str(xmax))
print("ymin = " + str(ymin))
print("ymax = " + str(ymax))
print("zmin = " + str(zmin))
print("zmax = " + str(zmax))
print("sim_velx_offset = ", sim_velx_offset)
print("sim_vely_offset = ", sim_vely_offset)
print("sim_velz_offset = ", sim_velz_offset)
'''
