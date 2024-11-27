import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
import glob
import my_flash_module as mym

sink_evol_pickle = sys.argv[1]
zoom_in_sink = sys.argv[2]
sim_dir = '/scratch/ek9/ccf100/sf_outflow/r1024mM5Ma2A1oh/'

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

form_time = yt.YTArray(sink_data[zoom_in_sink]['time'][0], 's')
form_position = yt.YTArray([sink_data[zoom_in_sink]['posx'][0], sink_data[zoom_in_sink]['posy'][1], sink_data[zoom_in_sink]['posy'][2]], 'cm')
box_length = yt.YTQuantity(0.1, 'pc')

xmin = form_position[0] - box_length/2
xmax = form_position[0] + box_length/2
ymin = form_position[1] - box_length/2
ymax = form_position[1] + box_length/2
zmin = form_position[2] - box_length/2
zmax = form_position[2] + box_length/2

#Find simfile right before sink forms
files = sorted(glob.glob(sim_dir + '*plt_cnt*'))
form_file = mym.find_files([form_time], files)
import pdb
pdb.set_trace()
