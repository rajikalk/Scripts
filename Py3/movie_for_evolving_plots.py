#!/usr/bin/env python
import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import os
from matplotlib import transforms
import glob
import my_module as mym
import pickle
import yt

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=5000, type=int)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def get_files(path, args):
    source_directory = sorted(glob.glob(path + 'WIND_' + slice + '*'))
    return source_directory

def has_sinks(file):
    try:
        part_file = file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        if f[list(f.keys())[1]][-1][-1] > 0:
            f.close()
            return True
        else:
            f.close()
            return False
    except:
        f = h5py.File(file, 'r')
        if "particlepositions" in list(f.keys()):
            f.close()
            return True
        else:
            f.close()
            return False

#=======MAIN=======
def main():
    # Read in directories:
    try:
        os.makedirs(save_dir)
    except:
        print("save directory exists")
    
    args = parse_inputs()
    mym.set_global_font_size(args.text_font)
    files = sorted(glob.glob(path + '*slice*'))

    #Define arrays
    time_array = []
    mass_array = []
    accr_array = []
    
    #get end time and set xlim
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=args.end_time)
    xlim = [0.0, args.end_time]
    xlabel = "Time ($yr$)"
    ylabel = "Accretion Rate ($M_\odot/yr$)"
    
    #Find sink formation time
    sink_formation_time = mym.find_sink_formation_time(files)
    
    frame_counter = 0
    for file in files:
        file_name = save_dir + "movie_frame_" + ("%06d" % frame_counter)
        frame_counter = frame_counter + 1
        f = h5py.File(file, 'rf')
        time = yt.YTQuantity(f['time'][0], 's').in_units('yr').value - sink_formation_time
        mass = yt.YTQuantity(np.sum(f['particlemasses'][:]), 'g').in_units('Msun').value
        if len(mass_array):
            accretion = np.nan
        else:
            accretion = (mass - mass_array[-1])/(time - time_array[-1])
        time_array.append(time)
        mass_array.append(mass)
        accr_array.append(accretion)
        
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
        ax.set_ylabel(yabel, labelpad=-20, fontsize=args.text_font)
        ax.set_xlim(xlim)
        ax.plot(time_array, accretion)
        plt.gca().set_aspect('equal')
        plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
        eps_image = Image.open(file_name + ".eps")
        eps_image.load(scale=4)
        eps_image.save(file_name + ".jpg")
        print('Created frame', frame_counter-1)

if __name__ == '__main__': main()


