import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
import csv

max_file = [2870, 2870]
it = 0
field = "dens"
type = ["proj", "slice"] #or "slice"
axis = ["xz", "xy"]
zoom = [False, True] #boolean
zoom_cell = [0,226]
annotate_vel_freq = [8,2]
start_frame = 143
source_directory = ["/Users/rajikak/Output/CircumbinaryOutFlow_0.25_lref_9/", "/Users/rajikak/Output/CircumbinaryOutFlow_0.25_lref_9/"]
directory = [source_directory[0] + "WIND_" + type[0] + "_" + axis[0] + "_000000", source_directory[1] + "WIND_" + type[1] + "_" + axis[1] + "_000000"]
year = 31557600
au = 1.496e13
m_times = []
its = [[],[]]
width_au = []

#load times:
with open('movie_times.txt', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        m_times.append(float(row[0]))
print "loaded movie frame times"

#load initial files:
X = [[],[]]
Y = [[],[]]
X_vel = [[],[]]
Y_vel = [[],[]]
for it in range(len(directory)):
    f = h5py.File(directory[it], 'r')
    xmin = f['minmax_xyz'][0][0]/au
    xmax = f['minmax_xyz'][0][1]/au
    width_au.append(xmax)
    if axis[it] == "xy":
        ymin = f['minmax_xyz'][1][0]/au
        ymax = f['minmax_xyz'][1][1]/au
    else:
        ymin = f['minmax_xyz'][2][0]/au
        ymax = f['minmax_xyz'][2][1]/au
    dim = np.shape(f[field + "_" + type[it] + "_" + axis[it]])[0]
    cl = (xmax-xmin)/(dim)
    x = np.arange(xmin+(zoom_cell[it]*cl), xmax-(zoom_cell[it])*cl, cl)
    y = np.arange(ymin+(zoom_cell[it]*cl), ymax-(zoom_cell[it])*cl, cl)
    x_vel = np.arange(xmin+(((annotate_vel_freq[it]/2.)+zoom_cell[it])*cl), xmax-(zoom_cell[it]*cl), cl*annotate_vel_freq[it])
    y_vel = np.arange(ymin+(((annotate_vel_freq[it]/2.)+zoom_cell[it])*cl), ymax-(zoom_cell[it]*cl), cl*annotate_vel_freq[it])
    X[it], Y[it] = np.meshgrid(x, y)
    X_vel[it], Y_vel[it] = np.meshgrid(x_vel, y_vel)
print "created meshs"

#find usable plots:
for file in range(len(its)):
    it = 0
    mit = 0
    diff_prev = np.inf
    while mit < len(m_times)-1:
        mtime = m_times[mit]
        while it < max_file[file]:
            directory = source_directory[file] + "WIND_" + type[file] + "_" + axis[file] + "_" + ("%06d" % it)
            f = h5py.File(directory, 'r')
            diff = abs(mtime - float(str(f['time'][0])))
            if diff == 0:
                its[file].append(it)
                mit = mit + 1
                if mit == len(m_times):
                    break
                else:
                    mtime = m_times[mit]
            elif np.isinf(diff_prev) and diff != 0:
                diff_prev = diff
            elif diff > diff_prev:
                its[file].append(it-1)
                diff_prev = np.inf
                it = it - 1
                mit = mit + 1
                if mit == len(m_times):
                    break
                else:
                    mtime = m_times[mit]
            elif diff < diff_prev:
                diff_prev = diff
            if it == max_file[file]-1 and mtime == m_times[-1]:
                its[file].append(it)
            it = it + 1
print "found usable movie plots"

#find extremes for colour bar
max = []
min = []
for file in range(len(its)):
    for it in its[file]:
        directory = source_directory[file] + "WIND_" + type[file] + "_" + axis[file] + "_" + ("%06d" % it)
        f = h5py.File(directory, 'r')
        if zoom[file]:
            for x in range(zoom_cell[file], dim-zoom_cell[file]):
                den_field = f[field + '_' + type[file] + '_' + axis[file]][x]
                if shape(den_field) == (dim,1):
                    den_field = den_field.transpose()[0]
                max_val = np.max(den_field[zoom_cell[file]:dim-zoom_cell[file]])
                min_val = np.min(den_field[zoom_cell[file]:dim-zoom_cell[file]])
                if type[file] == "proj":
                    if axis[file] == "xz":
                        max_val = max_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                        min_val = min_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                    else:
                        max_val = max_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
                        min_val = min_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
        else:
            max_val = np.max(f[field + '_' + type[file] + '_' + axis[file]])
            min_val = np.min(f[field + '_' + type[file] + '_' + axis[file]])
            if type[file] == "proj":
                if axis[file] == "xz":
                    max_val = max_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                    min_val = min_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                else:
                    max_val = max_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
                    min_val = min_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
        max.append(max_val)
        min.append(min_val)
    print "found maxes and mins for plot", file
max_str = '%.1e' % np.max(max)
min_str = '%.1e' % np.min(min)
max_round = str(0.5*floor(2.0 * float(max_str[0:3])))
min_round = str(0.5*ceil(2.0 * float(min_str[0:3])))
cbar_max = float(max_round+max_str[-4:])
cbar_min = float(min_round+min_str[-4:])
print "found colour bar extremes"

#plot figures
frame_val = start_frame
its = [its[0][start_frame:],its[1][start_frame:]]
for it in range(len(its[0])):
    velx = [[],[]]
    vely = [[],[]]
    image = [[],[]]
    magx = [[],[]]
    magy = [[],[]]
    part_pos= [[],[]]
    particles = [False, False]
    part_mass = [[],[]]
    for file in range(len(source_directory)):
        directory = source_directory[file] + "WIND_" + type[file] + "_" + axis[file] + "_" + ("%06d" % its[file][it])
        f = h5py.File(directory, 'r')
        for x in range(zoom_cell[file], dim-zoom_cell[file]):
            image_val = f[field + "_" + type[file] + "_" + axis[file]][x]
            magx_val = f["mag"+axis[file][0]+ "_" + type[file] + "_"+axis[file]][x]
            magy_val = f["mag"+axis[file][1]+ "_" + type[file] + "_"+axis[file]][x]
            if shape(image_val) == (dim,1):
                image_val = image_val.transpose()
                magx_val = magx_val.transpose()
                magy_val = magy_val.transpose()
            image_val = image_val[0][zoom_cell[file]: dim-zoom_cell[file]]
            if type[file] == "proj":
                if axis[file] == "xy":
                    image_val = image_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
                else:
                    image_val = image_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
            magx_val = magx_val[0][zoom_cell[file]: dim-zoom_cell[file]]
            magy_val = magy_val[0][zoom_cell[file]: dim-zoom_cell[file]]
            image[file].append(image_val)
            magx[file].append(magx_val)
            magy[file].append(magy_val)
        image[file] = np.array(image[file])
        magx[file] = np.array(magx[file])
        magy[file] = np.array(magy[file])
        for x in range(zoom_cell[file], dim-zoom_cell[file], annotate_vel_freq[file]):
            velx_temp = f["vel" + axis[file][0] + "_" + type[file] + "_" + axis[file]][x]
            vely_temp = f["vel" + axis[file][1] + "_" + type[file] + "_" + axis[file]][x]
            if shape(velx_temp) == (dim, 1):
                velx_temp = velx_temp.transpose()
                vely_temp = vely_temp.transpose()
            velx[file].append(velx_temp[0][zoom_cell[file]:dim-zoom_cell[file]:annotate_vel_freq[file]])
            vely[file].append(vely_temp[0][zoom_cell[file]:dim-zoom_cell[file]:annotate_vel_freq[file]])
        if len(f.keys()) > 12:
            particles[file] = True
            ypos_sort = np.sort(f["particlepositions"][1])
            ypos_sort[:] = ypos_sort[::-1]
            for ypos in ypos_sort:
                for pos in range(len(f["particlepositions"][1])):
                    if f["particlepositions"][1][pos] == ypos:
                        part_pos[file].append(f["particlepositions"][0]/au)
                        part_pos[file].append(f["particlepositions"][1]/au)
                        part_pos[file].append(f["particlepositions"][2]/au)
            part_mass_max = np.max(f["particlemasses"])
            part_mass[file] = f["particlemasses"]/part_mass_max
    plt.clf()
    fig, axes = plt.subplots(nrows=1,ncols=2)
    fig.subplots_adjust(wspace=0.0)
    fig.set_size_inches(14,6)
    axes[1].set_yticklabels(axes[1].get_yticklabels(), visible=False)
    for ax in range(len(axes)):
        plot = axes[ax].pcolormesh(X[ax], Y[ax], image[ax], cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max))
        Q = quiver(X_vel[ax], Y_vel[ax], velx[ax], vely[ax])
        axes[ax].streamplot(X[ax], Y[ax], magx[ax], magy[ax])
        if particles[ax]:
            if axis[ax] == "xz":
                axes[ax].scatter(part_pos[file][0], part_pos[file][2], c=part_mass[file], cmap=mpl.cm.gray)
            else:
                axes[ax].scatter(part_pos[file][0], part_pos[file][1], c=part_mass[file], cmap=mpl.cm.gray)
        axes[ax].set_xlabel('x (AU)')
        if zoom[ax]:
            axes[ax].set_xlim([ceil(X[ax][0][0]/10)*10, floor(X[ax][-1][-1]/10)*10])
            axes[ax].set_ylim([ceil(Y[ax][0][0]/10)*10, floor(Y[ax][-1][-1]/10)*10])
        else:
            axes[ax].set_xlim([-width_au[ax], width_au[ax]])
            axes[ax].set_ylim([-width_au[ax], width_au[ax]])
    cbar = fig.colorbar(plot, pad=0.0)
    cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    fig.suptitle('time='+str((f['time'][0])/year)+'years')
    if axis[0] == "xy":
        axes[0].set_ylabel('y (AU)', labelpad=-20)
    else:
        axes[0].set_ylabel('z (AU)', labelpad=-20)
    plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".png", bbox_inches='tight')
    frame_val = frame_val + 1
    print 'done frame_val =', frame_val