#!/usr/bin/env python
import sys
import yt
yt.enable_parallelism()
import glob
import my_flash_module as mym
import numpy as np
import pickle

input_dir = sys.argv[1]
save_dir = sys.argv[2]
plot_time = yt.YTQuantity(4220, 'yr')
pickle_file = save_dir + "time_" + str(int(plot_time.value)) +".pkl"
plot_width = 200
quiver_arrows = 1
axis = 'z'

x_image_min = yt.YTQuantity(-1*plot_width/2, 'au')
x_image_max = yt.YTQuantity(plot_width/2, 'au')
#x_image_min = -1*ds.domain_width.in_units('au')[0]/2
#x_image_max = ds.domain_width.in_units('au')[0]/2
x_range = np.linspace(x_image_min, x_image_max, 800)
X_image, Y_image = np.meshgrid(x_range, x_range)
annotate_space = (x_image_max - x_image_min)/quiver_arrows
x_ind = []
y_ind = []
counter = 0
while counter < quiver_arrows:
    val = annotate_space*counter + annotate_space/2. + x_image_min
    x_ind.append(int(val))
    y_ind.append(int(val))
    counter = counter + 1
X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)

files = sorted(glob.glob(input_dir + '*plt_cnt*'))
fn = mym.find_files([plot_time.value], files)[0]
part_file = 'part'.join(fn.split('plt_cnt'))
ds = yt.load(fn, particle_filename=part_file)
dd = ds.all_data()

#Altnerative calculation
Primary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[0], dd['particle_posy'].in_units('au')[0], dd['particle_posz'].in_units('au')[0]])
Secondary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[1], dd['particle_posy'].in_units('au')[1], dd['particle_posz'].in_units('au')[1]])
d_pos = Secondary_pos - Primary_pos

Primary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[0], dd['particle_vely'].in_units('km/s')[0], dd['particle_velz'].in_units('km/s')[0]])
Secondary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[1], dd['particle_vely'].in_units('km/s')[1], dd['particle_velz'].in_units('km/s')[1]])
d_vel = Secondary_vel - Primary_vel
L_vec = np.cross(d_pos, d_vel).T
L_vec_norm = L_vec/np.sqrt(np.sum(L_vec**2))

part_info = {'particle_mass':dd['particle_mass'].in_units('msun'),
             'particle_position':yt.YTArray([Primary_pos, Secondary_pos]).T,
             'particle_velocities':yt.YTArray([Primary_vel, Secondary_vel]).T,
             'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
             'particle_tag':dd['particle_tag'],
             'particle_form_time':dd['particle_creation_time']}

slice_field_list = [('flash', 'dens')]
slice_field_list = slice_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+axis not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+axis not in field[1])]

slice_dict = {}
for sto, field in yt.parallel_objects(slice_field_list, storage=slice_dict):
    #print("Projecting field", field, "on rank", rank)
    slc = yt.OffAxisSlicePlot(ds, L_vec_norm, field, width=plot_width, center=Primary_pos, north_vector=[0, 1, 0])
    slice_array = slc.frb.data[field].in_cgs()
    #print(field, "projection =", proj_array)
    sto.result_id = field[1]
    sto.result = slice_array

velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, slice_dict[list(slice_dict.keys())[1]], slice_dict[list(slice_dict.keys())[2]], no_of_quivers=quiver_arrows)
file = open(pickle_file, 'wb')
pickle.dump((X_image, Y_image, slice_dict[list(slice_dict.keys())[0]], slice_dict[list(slice_dict.keys())[3]], slice_dict[list(slice_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, plot_time), file)
file.close()
print("created pickle", pickle_file)

import pdb
pdb.set_trace()

