#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-22')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-19)
    parser.add_argument("-G_mass", "--Gas_mass", type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories
args = parse_inputs()

cbar_max = args.colourbar_max
try:
    cbar_min = float(args.colourbar_min)
except:
    cbar_min = float(args.colourbar_min[1:])

sys.stdout.flush()
CW.Barrier()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False and rank==0:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

sys.stdout.flush()
proj_pickles = sorted(glob.glob(save_dir + "movie_frame_*_proj.pkl"))
part_pickles = sorted(glob.glob(save_dir + "movie_frame_*_part.pkl"))

sys.stdout.flush()
CW.Barrier()

x = np.linspace(-2, 2, 800)
y = np.linspace(-2, 2, 800)
xlim = [x[0], x[-1]]
ylim = [y[0], y[-1]]
X, Y = np.meshgrid(x, y)

sys.stdout.flush()
CW.Barrier()

rit = -1
for pick_it in range(len(proj_pickles)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        save_name = save_dir + "movie_frame_" + ("%06d" % pick_it)
        if os.path.exists(save_name+".jpg") == False:
            file = open(proj_pickles[pick_it], 'rb')
            image, time_val = pickle.load(file)
            file.close()
            
            plt.clf()
            fig, ax = plt.subplots()
            ax.set_xlabel("$x$ (pc)", labelpad=-1, fontsize=args.text_font)
            ax.set_ylabel("$y$ (pc)", fontsize=args.text_font) #, labelpad=-20
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar(plot, pad=0.0)
            
            file = open(part_pickles[pick_it], 'rb')
            particle_x_pos, particle_y_pos, particle_mass = pickle.load(file)
            file.close()
            
            ax.scatter((particle_x_pos.value - 2), (particle_y_pos.value - 2), color='c', s=0.5)
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)

            SFE_val = np.round(100*(np.sum(particle_mass).value/args.Gas_mass), decimals=1)
            time_string = "$SFE$="+str(SFE_val)+"%"
            time_string_raw = r"{}".format(time_string)
            time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
            
            plt.savefig(save_name + ".jpg", format='jpg', bbox_inches='tight')
            plt.savefig(save_name + ".pdf", format='pdf', bbox_inches='tight')
            print('Created frame on rank', rank, 'at time of', str(time_val), 'to save_dir:', proj_pickles[pick_it] + '.jpg')
            
sys.stdout.flush()
CW.Barrier()

print("Finished making frames")
