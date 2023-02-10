#!/usr/bin/env python
import numpy as np
import sys
import pickle
import csv
import matplotlib.pyplot as plt

plot_pickles = ['/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_1_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_2_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/bound_core_frag_(106_77)_3_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_1_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_2_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/unbound_core_frag_121_104_3_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_1_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_2_part.pkl', '/lustre/astro/rlk/Movie_frames/Ramses/Global/G100/256/XY/Formation_pathways/dynamical_capt_(101_(13_[77_106]))_3_part.pkl']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(two_col_width,two_col_width))

for pick_it in range(len(plot_pickles)):
    file = open(pickle_file, 'rb')
    particle_x_pos, particle_y_pos, particle_masses, max_seps, sink_creation_time_pick, center_pos = pickle.load(file)
    file.close()
    
    file = open("".join(pickle_file.split('_part')), 'rb')
    image, time_val = pickle.load(file)
    file.close()
    
import pdb
pdb.set_trace()

