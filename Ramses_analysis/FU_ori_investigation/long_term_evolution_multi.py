import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt
import sys


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_id", help="which sink?", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}" "\sisetup{detect-all}" r"\usepackage{helvet}" r"\usepackage{sansmath}" "\sansmath"               # <- tricky! -- gotta actually tell tex to use!

sink_inds = [17, 45, 51, 71, 75, 85, 101, 103, 176, 177, 258, 272, 292]

for sink_ind in sink_inds:
    pickle_file = '/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_'+str(sink_ind)+'.pkl'
    if os.path.isfile(pickle_file):
        file_open = open(pickle_file, 'rb')
        particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
        file_open.close()
        
    mass_ratio = 
