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
proj_thickness = 2
colourbar_min = None
colourbar_max = None
standard_vel = None

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


part_info = mym.get_particle_data(ds, sink_id=sink_id)
pos_array = yt.YTArray([Primary_pos, Secondary_pos])
east_unit_vector = [1, 0, 0]
north_unit = [0, 1, 0]

projected_particle_posy = projected_vector(pos_array, north_unit)
proj_part_y_mag = np.sqrt(np.sum((projected_particle_posy**2), axis=1))
proj_part_y_unit = (projected_particle_posy.T/proj_part_y_mag).T
north_sign = np.dot(north_unit, proj_part_y_unit.T)
proj_part_y = proj_part_y_mag*north_sign
proj_part_y = np.nan_to_num(proj_part_y)

projected_particle_posx = projected_vector(pos_array, east_unit_vector)
proj_part_x_mag = np.sqrt(np.sum((projected_particle_posx**2), axis=1))
proj_part_x_unit = (projected_particle_posx.T/proj_part_x_mag).T
east_sign = np.dot(east_unit_vector, proj_part_x_unit.T)
proj_part_x = proj_part_x_mag*east_sign
proj_part_x = np.nan_to_num(proj_part_x)

part_info['particle_position'] = np.array([[proj_part_x[0].value, proj_part_x[1].value],[proj_part_y[0].value, proj_part_y[1].value]])

'''
part_info = {'particle_mass':dd['particle_mass'].in_units('msun'),
             'particle_position':yt.YTArray([Primary_pos, Secondary_pos]).T,
             'particle_velocities':yt.YTArray([Primary_vel, Secondary_vel]).T,
             'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
             'particle_tag':dd['particle_tag'],
             'particle_form_time':dd['particle_creation_time']}
'''

proj_field_list = [('flash', 'dens')]
proj_field_list = proj_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+axis not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+axis not in field[1])]

proj_dict = {}
for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
    #print("Projecting field", field, "on rank", rank)
    proj = yt.OffAxisProjectionPlot(ds, L_vec_norm, field, width=plot_width, center=Primary_pos, north_vector=[0, 1, 0], depth=(proj_thickness, 'AU'))
    proj_array = proj.frb.data[field].in_cgs()/yt.YTQuantity(proj_thickness, 'AU').in_units('cgs')
    #print(field, "projection =", proj_array)
    sto.result_id = field[1]
    sto.result = proj_array

velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=quiver_arrows)
file = open(pickle_file, 'wb')
pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, plot_time), file)
file.close()
print("created pickle", pickle_file)

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_flash_module as mym

file = open(pickle_file, 'rb')
X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
file.close()

primary_ind = np.argmin(part_info['particle_form_time'])
part_velx = part_info['particle_velocities'][0][primary_ind]
part_vely = part_info['particle_velocities'][1][primary_ind]
velx = velx - part_velx.value
vely = vely - part_vely.value

plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('AU', labelpad=-1, fontsize=10)
ax.set_ylabel('AU', fontsize=10) #, labelpad=-20
xlim = [np.min(X_image).value, np.max(X_image).value]
ylim = [np.min(Y_image).value, np.max(Y_image).value]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

if colourbar_min == None:
    cmin = np.min(image)
else:
    cmin = colourbar_min
    
if colourbar_max == None:
    cmax = np.max(image)
else:
    cmax = colourbar_max
cbar_lims = [cmin, cmax]

if standard_vel == None:
    stdvel = 5
else:
    stdvel = standard_vel

cmap=plt.cm.gist_heat
plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
plt.gca().set_aspect('equal')

plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
cbar = plt.colorbar(plot, pad=0.0)

try:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
except:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)

primary_ind = np.argmax(part_info['particle_mass'])
part_info['particle_position'][0] = part_info['particle_position'][0] - part_info['particle_position'][0][primary_ind]
part_info['particle_position'][1] = part_info['particle_position'][1] - part_info['particle_position'][1][primary_ind]                    #Update particle position
mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=7)

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')
    
cbar.set_label("Density (" + str(image.units)+")", rotation=270, labelpad=14, size=10)

time_string = "$t$="+str(int(time_val))+"yr"
time_string_raw = r"{}".format(time_string)
time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

file_name = "Time_" + str(int(plot_time.value)) +".jpg"
plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
