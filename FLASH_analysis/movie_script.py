#!/usr/bin/env python
import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import pickle

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
X_vel, Y_vel = np.meshgrid(x_ind, y_ind)

#Now let's iterate over the files and get the images we want to plot
file_counter = -1
for fn in yt.parallel_objects(movie_files, njobs=int(size/5)):
    file_counter = file_counter + 1
    pickle_file = output_dir+"movie_frame_"+("%06d" % file_counter)+".pkl"
    proj_root_rank = int(rank/5)*5
    part_file = 'part'.join(fn.split('plt_cnt'))
    ds = yt.load(fn, particle_filename=part_file)
    
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
                 'particle_position':positions}
    else:
        has_particles = False
        part_info = {}
        
    file = open(pickle_file, 'wb')
    pickle.dump((X_image, Y_image, proj_dict[proj_field_list[0][1]], proj_dict[proj_field_list[3][1]], proj_dict[proj_field_list[4][1]], X_image_vel, Y_image_vel, proj_dict[proj_field_list[1][1]], proj_dict[proj_field_list[2][1]], part_info), file)
    file.close()
