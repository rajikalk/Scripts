import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
from pyramses import rsink
import sys
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
#================================================================================
args = parse_inputs()

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)
    
#Define empty particle array:
particle_data = {}
particle_data.update({'particle_tag':[]})
particle_data.update({'time':[]})
particle_data.update({'posx':[]})
particle_data.update({'posy':[]})
particle_data.update({'posz':[]})
particle_data.update({'velx':[]})
particle_data.update({'vely':[]})
particle_data.update({'velz':[]})
particle_data.update({'anglx':[]})
particle_data.update({'angly':[]})
particle_data.update({'anglz':[]})
particle_data.update({'mass':[]})
particle_data.update({'mdot':[]})
sink_form_time = 0
    
loaded_sink_data = rsink(datadir=path, all=True)

