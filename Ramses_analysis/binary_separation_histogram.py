import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
from pyramses import rsink
import sys
import os
import yt

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

ul = 4.*3600.*180./np.pi
lb = 361
ub = 1724
is1 = args.sink_number
is2 = is1 + 1

ss = []
dist = []
for io in range(lb,ub):
  s=rsink(io,datadir=path)
  d=((s['x'][is1] - s['x'][is2])**2 +
     (s['y'][is1] - s['y'][is2])**2 +
     (s['z'][is1] - s['z'][is2])**2)**0.5*ul
  dist.append(d)
  ss.append(s)
  
file = open('separation_data.pkl', 'wb')
pickle.dump((dist, ss), file)
file.close()

plt.figure()
plt.hist(dist)
plt.savefig('hist.pdf')

