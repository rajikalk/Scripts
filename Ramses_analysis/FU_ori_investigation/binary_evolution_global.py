import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt
from matplotlib.ticker import FormatStrFormatter

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

time_window = yt.YTQuantity(100, 'yr')

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

age = yt.YTArray(particle_data['time']).in_units('kyr')
mass_ratio = yt.YTArray(particle_data['mass']).T[1]/yt.YTArray(particle_data['mass']).T[0]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(two_col_width, 1.5*two_col_width), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.1)

axs.flatten()[0].plot(age, yt.YTArray(particle_data['mass']).T[0], label='Primary')
axs.flatten()[0].plot(age, yt.YTArray(particle_data['mass']).T[1], label='Secondary')
axs.flatten()[0].set_xlim()

axs.flatten()[1].plot(age, mass_ratio)

axs.flatten()[2].plot(age, particle_data['separation'])

axs.flatten()[3].plot(age, particle_data['eccentricity'])


