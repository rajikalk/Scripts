#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import pickle
import os

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="z")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
args = parse_inputs()

#Get movie files
movie_files = sorted(glob.glob(input_dir + '*_plt_cnt*'))

#Calculate image grid:
fn = movie_files[-1]
part_file = 'part'.join(fn.split('plt_cnt'))
ds = yt.load(fn, particle_filename=part_file)
x_image_min = -1*ds.domain_width.in_units('au')[0]/2
x_image_max = ds.domain_width.in_units('au')[0]/2
x_range = np.linspace(x_image_min, x_image_max, 800)
X_image, Y_image = np.meshgrid(x_range, x_range)
annotate_space = (x_image_min - x_image_max)/32.
x_ind = []
y_ind = []
counter = 0
while counter < 32:
    val = annotate_space*counter + annotate_space/2. + x_image_min
    x_ind.append(int(val))
    y_ind.append(int(val))
    counter = counter + 1
X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)

#Now let's iterate over the files and get the images we want to plot
file_counter = -1
for fn in yt.parallel_objects(movie_files, njobs=int(size/5)):
    file_counter = file_counter + 1
    pickle_file = output_dir+"movie_frame_"+("%06d" % file_counter)+".pkl"
    make_pickle = False
    if os.path.isfile(pickle_file) == False:
        make_pickle = True
    if make_pickle:
        proj_root_rank = int(rank/5)*5
        part_file = 'part'.join(fn.split('plt_cnt'))
        ds = yt.load(fn, particle_filename=part_file)
        time_val = ds.current_time.in_units('yr')
        
        #make list of projection fields: density, velocity, magnetic field
        proj_field_list = [('flash', 'dens')] + \
            [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+args.axis not in field[1])] + \
            [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+args.axis not in field[1])]
        
        #This is the dictionary where the projected arrays will be saved:
        proj_dict = {}
        for field in proj_field_list:
            proj_dict.update({field[1]:[]})
        
        #Make projections of each field
        for field in yt.parallel_objects(proj_field_list):
            proj = yt.ProjectionPlot(ds, args.axis, field, method='integrate')
            thickness = (proj.bounds[1] - proj.bounds[0]).in_cgs() #MIGHT HAVE TO UPDATE THIS LATER
            proj_array = proj.frb.data[field].in_cgs()/thickness
            if rank == proj_root_rank:
                proj_dict[field[1]] = proj_array
            else:
                file = open(pickle_file.split('.pkl')[0] + '_proj_data_' + str(proj_root_rank)+ str(proj_dict_keys.index(field[1])) + '.pkl', 'wb')
                pickle.dump((field[1], proj_array), file)
                file.close()
            
        #gather projection arrays
        if rank == proj_root_rank and size > 1:
            for kit in range(1,len(proj_dict_keys)):
                file = open(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl', 'rb')
                key, proj_array = pickle.load(file)
                file.close()
                proj_dict[key] = proj_array
                os.remove(pickle_file.split('.pkl')[0] + '_proj_data_' +str(proj_root_rank) +str(kit)+'.pkl')

        #Get particle data:
        dd = ds.all_data()
        if len([field for field in ds.field_list if 'particle_mass' in field[1]]) > 0:
            has_particles = True
            part_mass = dd['particle_mass'].in_units('msun')
            part_pos_fields = [field for field in ds.field_list if ('particle_pos' in field[1])&(field[0]=='all')&(field[1]!='particle_pos'+args.axis)]
            part_pos_x = dd[part_pos_fields[0]].in_units('au')
            part_pos_y = dd[part_pos_fields[1]].in_units('au')
            positions = np.array([part_pos_x,part_pos_y])
            part_info = {'particle_mass':part_mass,
                     'particle_position':positions,
                     'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
                     'particle_tag':dd['particle_tag']}
        else:
            has_particles = False
            part_info = {}
            
        file = open(pickle_file, 'wb')
        pickle.dump((X_image, Y_image, proj_dict[proj_field_list[0][1]], proj_dict[proj_field_list[3][1]], proj_dict[proj_field_list[4][1]], X_image_vel, Y_image_vel, proj_dict[proj_field_list[1][1]], proj_dict[proj_field_list[2][1]], part_info, time_val), file)
        file.close()

#Make frames.
import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_flash_module as mym

#Let's get the pickle files
pickle_files = sorted(glob.glob(output_dir+"movie_frame_*.pkl"))

rit = -1
for pickle_file in pickle_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        print('making frame from', pickle_file, 'on rank', rank)
        frame_no = int(pickle_file.split('_')[-1].split('.')[0])
        file = open(pickle_file, 'rb')
        X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
        file.close()

        file_name = output_dir + "movie_frame_" + ("%06d" % frame_no)
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_xlabel('AU', labelpad=-1, fontsize=10)
        ax.set_ylabel('AU', fontsize=10) #, labelpad=-20
        xlim = [np.min(X_image), np.max(X_image)]
        ylim = [np.min(Y_image), np.max(Y_image)]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        cmap=plt.cm.gist_heat
        plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
        plt.gca().set_aspect('equal')

        if frame_no > 0 or time_val > -1.0:
            plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
        else:
            plt.streamplot(X_image, Y_image, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
        cbar = plt.colorbar(plot, pad=0.0)
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], standard_vel=5, Z_val=None)

        if len(part_info.keys())>0:
            mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7)

        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
        plt.tick_params(axis='both', which='major')# labelsize=16)
        for line in ax.xaxis.get_ticklines():
            line.set_color('white')
        for line in ax.yaxis.get_ticklines():
            line.set_color('white')
        if args.debug_plotting != 'False':
            plt.savefig("Test_826.jpg", format='jpg', bbox_inches='tight')

        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

        if size > 1:
            try:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
                plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
            except:
                print("couldn't save for the dviread.py problem. Make frame " + str(frame_no) + " on ipython")
        else:
            plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
            plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
            print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
