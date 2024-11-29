import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
import glob
import my_flash_module as mym
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_flash_module as mym
import matplotlib.patches as patches

sink_evol_pickle = sys.argv[1]
primary_sink = sys.argv[2]
secondary_sink = sys.argv[3]
sim_dir = '/scratch/ek9/ccf100/sf_outflow/r1024mM5Ma2A1oh/'

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-box", "--box_length", help="how big do you want your resimulate box?", type=float, default=0.1)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()
box_length = yt.YTQuantity(args.box_length, 'pc')

file = open(sink_evol_pickle, 'rb')
sink_data, prev_line_counter = pickle.load(file)
file.close()

#Get binary CoM when secondary have 0.2Msun

secondary_mass_ind = np.argmin(abs(yt.YTArray(sink_data[secondary_sink]['mass'], 'g').in_units('Msun').value - 0.2))
binary_characteristic_time = yt.YTQuantity(sink_data[secondary_sink]['time'][secondary_mass_ind], 's')
binary_masses = yt.YTArray([sink_data[primary_sink]['mass'][secondary_mass_ind], sink_data[secondary_sink]['mass'][secondary_mass_ind]], 'g')
binary_positions = yt.YTArray([[sink_data[primary_sink]['posx'][secondary_mass_ind], sink_data[secondary_sink]['posx'][secondary_mass_ind]], [sink_data[primary_sink]['posy'][secondary_mass_ind], sink_data[secondary_sink]['posy'][secondary_mass_ind]], [sink_data[primary_sink]['posz'][secondary_mass_ind], sink_data[secondary_sink]['posz'][secondary_mass_ind]]], 'cm')
binary_velocities = yt.YTArray([[sink_data[primary_sink]['velx'][secondary_mass_ind], sink_data[secondary_sink]['velx'][secondary_mass_ind]], [sink_data[primary_sink]['vely'][secondary_mass_ind], sink_data[secondary_sink]['vely'][secondary_mass_ind]], [sink_data[primary_sink]['velz'][secondary_mass_ind], sink_data[secondary_sink]['velz'][secondary_mass_ind]]], 'cm/s')
CoM_pos = (binary_positions.T[0] * binary_masses[0] + binary_positions.T[1] * binary_masses[1])/np.sum(binary_masses)
CoM_vel = (binary_velocities.T[0] * binary_masses[0] + binary_velocities.T[1] * binary_masses[1])/np.sum(binary_masses)


#Find frame before primary formation
form_time = yt.YTArray(sink_data[primary_sink]['time'][0], 's')
files = sorted(glob.glob(sim_dir + '*plt_cnt*'))
files = [ x for x in files if "_proj_" not in x ]
form_file = mym.find_files([form_time.value], files)
prev_file = files[files.index(form_file[0])-1]

#Calculate bulk file
part_file = 'part'.join(prev_file.split('plt_cnt'))
ds = yt.load(prev_file, particle_filename=part_file)

dt = binary_characteristic_time.in_units('yr') - ds.current_time.in_units('yr')
shifted_CoM = (CoM_pos + ((-1*CoM_vel)*dt.in_units('s'))).in_units('au')
shifted_CoM = yt.YTArray([sink_data[primary_sink]['posx'][0], sink_data[primary_sink]['posy'][0], sink_data[primary_sink]['posz'][0]], 'cm').in_units('au')

xmin = shifted_CoM[0] - (box_length.in_units('au')/2)
xmax = shifted_CoM[0] + (box_length.in_units('au')/2)
ymin = shifted_CoM[1] - (box_length.in_units('au')/2)
ymax = shifted_CoM[1] + (box_length.in_units('au')/2)
zmin = shifted_CoM[2] - (box_length.in_units('au')/2)
zmax = shifted_CoM[2] + (box_length.in_units('au')/2)

#Calculate bulk velocity
left_corner = yt.YTArray([xmin, ymin, zmin], 'au')
right_corner = yt.YTArray([xmax, ymax, zmax], 'au')
region = ds.box(left_corner, right_corner)

sim_velx_offset = -1 * np.mean(region['velx'])
sim_vely_offset = -1 * np.mean(region['vely'])
sim_velz_offset = -1 * np.mean(region['velz'])

#Print results
print("sim_input_file = ", prev_file)
print("xmin = " + str(xmin))
print("xmax = " + str(xmax))
print("ymin = " + str(ymin))
print("ymax = " + str(ymax))
print("zmin = " + str(zmin))
print("zmax = " + str(zmax))
print("sim_velx_offset = ", sim_velx_offset)
print("sim_vely_offset = ", sim_vely_offset)
print("sim_velz_offset = ", sim_velz_offset)

#Zoom plots
x_range = np.linspace(xmin.in_units('pc'), xmax.in_units('pc'), 800)
y_range = np.linspace(ymin.in_units('pc'), ymax.in_units('pc'), 800)
X_image, Y_image = np.meshgrid(x_range, y_range)
annotate_space = (x_range[-1] - x_range[0])/32
x_ind = []
y_ind = []
counter = 0
while counter < 32:
    xval = annotate_space*counter + annotate_space/2. + x_range[0]
    yval = annotate_space*counter + annotate_space/2. + y_range[0]
    x_ind.append(float(xval))
    y_ind.append(float(val))
    counter = counter + 1
X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)
time_val = ds.current_time.in_units('yr').value





#Make projections
if os.path.exists('xy_proj_zoom.pkl') == False:
    x_range = np.linspace(ds.domain_left_edge[0].in_units('pc'), ds.domain_right_edge[0].in_units('pc'), 800)
    X_image, Y_image = np.meshgrid(x_range, x_range)
    x_image_min = ds.domain_left_edge[0].in_units('pc')
    x_image_max = ds.domain_right_edge[0].in_units('pc')
    annotate_space = (x_image_max - x_image_min)/32
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 32:
        val = annotate_space*counter + annotate_space/2. + x_image_min
        x_ind.append(float(val))
        y_ind.append(float(val))
        counter = counter + 1
    X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)
    time_val = ds.current_time.in_units('yr').value

    proj_field_list = [('flash', 'dens')]
    proj_field_list = proj_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('velz' not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('magz' not in field[1])]

    part_info = {'particle_mass':region['particle_mass'].in_units('msun'),
                 'particle_position':yt.YTArray([region['particle_posx'].in_units('pc'),region['particle_posy'].in_units('pc')]),
                 'particle_velocities':yt.YTArray([region['particle_velx'].in_units('cm/s'),region['particle_vely'].in_units('cm/s')]),
                 'accretion_rad':2.5*np.min(region['dx'].in_units('pc')),
                 'particle_tag':region['particle_tag'],
                 'particle_form_time':region['particle_creation_time']}
                 
    proj_dict = {}
    for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
        proj = yt.ProjectionPlot(ds, "z", field, method='integrate', data_source=region)
        proj_array = proj.frb.data[field].in_cgs()/box_length.in_units('cm')
        sto.result_id = field[1]
        sto.result = proj_array
    velx, vely, velz = mym.get_quiver_arrays(x_image_min.value, x_image_max.value, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=32)

    file = open('xy_proj_zoom.pkl', 'wb')
    pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
    file.close()

file = open('xy_proj_zoom.pkl', 'rb')
X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
file.close()

plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('x (AU)', labelpad=-1, fontsize=10)
ax.set_ylabel('y (AU)', fontsize=10) #, labelpad=-20
xlim = [np.min(X_image).value, np.max(X_image).value]
ylim = [np.min(Y_image).value, np.max(Y_image).value]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

cmin = np.min(image)
cmax = np.max(image)
cbar_lims = [cmin, cmax]
stdvel = 2
cmap=plt.cm.gist_heat
plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
plt.gca().set_aspect('equal')
plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, minlength=0.5, zorder=2)
cbar = plt.colorbar(plot, pad=0.0)
try:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
except:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
ax.scatter(part_info['particle_position'][0], part_info['particle_position'][1], color='c', s=0.5)

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')
    
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
time_string = "$t$="+str(int(time_val))+"yr"
time_string_raw = r"{}".format(time_string)
time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

square = patches.Rectangle((xmin.in_units('pc'), ymin.in_units('pc')), box_length, box_length, edgecolor='green', facecolor='none')
ax.add_patch(square)
CS = ax.contour(X_image,X_image,image, locator=plt.LogLocator(), linewidths=0.5, colors='blue', levels=[1.e-18])

plt.savefig("xy_proj_zoom.jpg", format='jpg', bbox_inches='tight', dpi=300)




















#For xy projection centred on z=shifted_CoM[2]
if os.path.exists('xy_proj.pkl') == False:
    proj_field_list = [('flash', 'dens')]
    proj_field_list = proj_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('velz' not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('magz' not in field[1])]

    proj_left_corner = yt.YTArray([ds.domain_left_edge[0], ds.domain_left_edge[1], shifted_CoM[2].in_units('cm')-box_length.in_units('cm')/2], 'cm')
    proj_right_corner = yt.YTArray([ds.domain_right_edge[0], ds.domain_right_edge[1], shifted_CoM[2].in_units('cm')+box_length.in_units('cm')/2], 'cm')
    proj_region = ds.box(proj_left_corner, proj_right_corner)
    part_info = {'particle_mass':proj_region['particle_mass'].in_units('msun'),
                 'particle_position':yt.YTArray([proj_region['particle_posx'].in_units('pc'),proj_region['particle_posy'].in_units('pc')]),
                 'particle_velocities':yt.YTArray([proj_region['particle_velx'].in_units('cm/s'),proj_region['particle_vely'].in_units('cm/s')]),
                 'accretion_rad':2.5*np.min(proj_region['dx'].in_units('pc')),
                 'particle_tag':proj_region['particle_tag'],
                 'particle_form_time':proj_region['particle_creation_time']}
                 
    proj_dict = {}
    for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
        proj = yt.ProjectionPlot(ds, "z", field, method='integrate', data_source=proj_region)
        proj_array = proj.frb.data[field].in_cgs()/box_length.in_units('cm')
        sto.result_id = field[1]
        sto.result = proj_array
    velx, vely, velz = mym.get_quiver_arrays(x_image_min.value, x_image_max.value, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=32)

    file = open('xy_proj.pkl', 'wb')
    pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
    file.close()

file = open('xy_proj.pkl', 'rb')
X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
file.close()

plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('x (AU)', labelpad=-1, fontsize=10)
ax.set_ylabel('y (AU)', fontsize=10) #, labelpad=-20
xlim = [np.min(X_image).value, np.max(X_image).value]
ylim = [np.min(Y_image).value, np.max(Y_image).value]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

cmin = 1.e-23
cmax = 1.e-17
cbar_lims = [cmin, cmax]
stdvel = 2
cmap=plt.cm.gist_heat
plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
plt.gca().set_aspect('equal')
plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, minlength=0.5, zorder=2)
cbar = plt.colorbar(plot, pad=0.0)
try:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
except:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
ax.scatter(part_info['particle_position'][0], part_info['particle_position'][1], color='c', s=0.5)

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')
    
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
time_string = "$t$="+str(int(time_val))+"yr"
time_string_raw = r"{}".format(time_string)
time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

square = patches.Rectangle((xmin.in_units('pc'), ymin.in_units('pc')), box_length, box_length, edgecolor='green', facecolor='none')
ax.add_patch(square)
CS = ax.contour(X_image,X_image,image, locator=plt.LogLocator(), linewidths=0.5, colors='blue', levels=[1.e-18])

plt.savefig("xy_proj.jpg", format='jpg', bbox_inches='tight', dpi=300)


#For xz projection centred on y=shifted_CoM[1]
if os.path.exists('xz_proj.pkl') == False:
    proj_field_list = [('flash', 'dens')]
    proj_field_list = proj_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vely' not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('magy' not in field[1])]

    proj_left_corner = yt.YTArray([ds.domain_left_edge[0], shifted_CoM[1].in_units('cm')-box_length.in_units('cm')/2, ds.domain_left_edge[2]], 'cm')
    proj_right_corner = yt.YTArray([ds.domain_right_edge[0], shifted_CoM[1].in_units('cm')+box_length.in_units('cm')/2, ds.domain_right_edge[2]], 'cm')
    proj_region = ds.box(proj_left_corner, proj_right_corner)
    part_info = {'particle_mass':proj_region['particle_mass'].in_units('msun'),
                 'particle_position':yt.YTArray([proj_region['particle_posx'].in_units('pc'),proj_region['particle_posz'].in_units('pc')]),
                 'particle_velocities':yt.YTArray([proj_region['particle_velx'].in_units('cm/s'),proj_region['particle_velz'].in_units('cm/s')]),
                 'accretion_rad':2.5*np.min(proj_region['dx'].in_units('pc')),
                 'particle_tag':proj_region['particle_tag'],
                 'particle_form_time':proj_region['particle_creation_time']}
                 
    proj_dict = {}
    for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
        proj = yt.ProjectionPlot(ds, "y", field, method='integrate', data_source=proj_region)
        proj_array = proj.frb.data[field].in_cgs()/box_length.in_units('cm')
        proj_array = proj_array.T
        sto.result_id = field[1]
        sto.result = proj_array
    velx, vely, velz = mym.get_quiver_arrays(x_image_min.value, x_image_max.value, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=32)

    file = open('xz_proj.pkl', 'wb')
    pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
    file.close()

file = open('xz_proj.pkl', 'rb')
X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
file.close()

plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('x (AU)', labelpad=-1, fontsize=10)
ax.set_ylabel('z (AU)', fontsize=10) #, labelpad=-20
xlim = [np.min(X_image).value, np.max(X_image).value]
ylim = [np.min(Y_image).value, np.max(Y_image).value]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

cmap=plt.cm.gist_heat
plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
plt.gca().set_aspect('equal')
plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, minlength=0.5, zorder=2)
cbar = plt.colorbar(plot, pad=0.0)
try:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
except:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
ax.scatter(part_info['particle_position'][0], part_info['particle_position'][1], color='c', s=0.5)

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')
    
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
time_string = "$t$="+str(int(time_val))+"yr"
time_string_raw = r"{}".format(time_string)
time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

square = patches.Rectangle((xmin.in_units('pc'), zmin.in_units('pc')), box_length, box_length, edgecolor='green', facecolor='none')
ax.add_patch(square)
#plt.contour([X_image.value, Y_image.value], image, [1.e-18])

plt.savefig("xz_proj.jpg", format='jpg', bbox_inches='tight', dpi=300)

#For yz projection centred on y=shifted_CoM[0]
if os.path.exists('yz_proj.pkl') == False:
    proj_field_list = [('flash', 'dens')]
    proj_field_list = proj_field_list + [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('velx' not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('magx' not in field[1])]

    proj_left_corner = yt.YTArray([shifted_CoM[0].in_units('cm')-box_length.in_units('cm')/2, ds.domain_left_edge[1], ds.domain_left_edge[2]], 'cm')
    proj_right_corner = yt.YTArray([shifted_CoM[0].in_units('cm')+box_length.in_units('cm')/2, ds.domain_right_edge[1], ds.domain_right_edge[2]], 'cm')
    proj_region = ds.box(proj_left_corner, proj_right_corner)
    part_info = {'particle_mass':proj_region['particle_mass'].in_units('msun'),
                 'particle_position':yt.YTArray([proj_region['particle_posy'].in_units('pc'),proj_region['particle_posz'].in_units('pc')]),
                 'particle_velocities':yt.YTArray([proj_region['particle_vely'].in_units('cm/s'),proj_region['particle_velz'].in_units('cm/s')]),
                 'accretion_rad':2.5*np.min(proj_region['dx'].in_units('pc')),
                 'particle_tag':proj_region['particle_tag'],
                 'particle_form_time':proj_region['particle_creation_time']}
                 
    proj_dict = {}
    for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
        proj = yt.ProjectionPlot(ds, "x", field, method='integrate', data_source=proj_region)
        proj_array = proj.frb.data[field].in_cgs()/box_length.in_units('cm')
        sto.result_id = field[1]
        sto.result = proj_array
    velx, vely, velz = mym.get_quiver_arrays(x_image_min.value, x_image_max.value, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=32)

    file = open('yz_proj.pkl', 'wb')
    pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
    file.close()

file = open('yz_proj.pkl', 'rb')
X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
file.close()

plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('y (AU)', labelpad=-1, fontsize=10)
ax.set_ylabel('z (AU)', fontsize=10) #, labelpad=-20
xlim = [np.min(X_image).value, np.max(X_image).value]
ylim = [np.min(Y_image).value, np.max(Y_image).value]
ax.set_xlim(xlim)
ax.set_ylim(ylim)

cmap=plt.cm.gist_heat
plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
plt.gca().set_aspect('equal')
plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, minlength=0.5, zorder=2)
cbar = plt.colorbar(plot, pad=0.0)
try:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
except:
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
ax.scatter(part_info['particle_position'][0], part_info['particle_position'][1], color='c', s=0.5)

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')
    
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
time_string = "$t$="+str(int(time_val))+"yr"
time_string_raw = r"{}".format(time_string)
time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

square = patches.Rectangle((ymin.in_units('pc'), zmin.in_units('pc')), box_length, box_length, edgecolor='green', facecolor='none')
ax.add_patch(square)
#plt.contour([X_image.value, Y_image.value], image, [1.e-18])

plt.savefig("yz_proj.jpg", format='jpg', bbox_inches='tight', dpi=300)
