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
from PIL import Image

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="mdot")
    parser.add_argument("-y_times", "--y_multiplier", help="What field to you wish to plot?", default=10000)
    parser.add_argument("-ylab", "--y_label", help="What is the label you want for the yaxis??", default="Accretion Rate ($10^{-4}M_\odot/yr$)")
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
    pickle_file = sys.argv[1]
    save_dir = sys.argv[2]
    try:
        os.makedirs(save_dir)
    except:
        print("save directory exists")
    
    args = parse_inputs()
    mym.set_global_font_size(args.text_font)
    #files = sorted(glob.glob(path + '*slice*'))

    
    #get end time and set xlim
    m_times = mym.generate_frame_times(pickle_file, args.time_step, presink_frames=0, end_time=args.end_time)
    m_times = m_times[args.start_frame:]
    #usable_files = mym.find_files(m_times, files)
    xlim = [0.0, args.end_time]
    xlabel = "Time ($yr$)"
    ylabel = args.y_label
    
    #Find sink formation time
    #sink_formation_time = mym.find_sink_formation_time(files)
    
    #Open pickle
    file_open = open(pickle_file, 'rb')
    particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
    file_open.close()
    
    #Smooth accretion
    window = 10 #in units years
    smoothed_time = []
    smoothed_quantity = []
    particle_data_time = np.array(particle_data['time'])
    particle_data_quantity = np.array(particle_data[args.field])
    for ind in range(len(particle_data['mdot'])):
        smoothing_inds = np.where((particle_data_time > particle_data['time'][ind]-window/2.)&(particle_data_time < particle_data['time'][ind]+window/2.))[0]
        time = np.mean(particle_data_time[smoothing_inds])
        quantity = np.mean(particle_data_quantity[smoothing_inds])
        smoothed_time.append(time)
        smoothed_quantity.append(quantity)
        
    for m_time in m_times:
        file_name = save_dir + "movie_frame_" + ("%06d" % (args.start_frame + m_times.index(m_time)))
        
        plot_ind = np.argmin(abs(np.array(smoothed_time) - m_time))
        
        plt.clf()
        plt.plot(smoothed_time, np.array(smoothed_quantity)*args.y_multiplier)
        plt.plot(smoothed_time[plot_ind], np.array(smoothed_quantity[plot_ind])*args.y_multiplier, 'ro')
        #plt.ticklabel_format(axis='y', style='sci')
        plt.xlabel(xlabel, labelpad=-1, fontsize=args.text_font)
        plt.ylabel(ylabel, labelpad=-2, fontsize=args.text_font)
        plt.ylim(bottom=0)
        plt.tight_layout()
        plt.savefig(file_name + ".eps", format='eps')#, bbox_inches='tight')#, pad_inches = 0.02)
        eps_image = Image.open(file_name + ".eps")
        eps_image.load(scale=4)
        eps_image.save(file_name + ".jpg")
        print('Created frame:', file_name)
    
if __name__ == '__main__': main()


