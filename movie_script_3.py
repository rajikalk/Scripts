import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

max_file = 1394
it = 0
field = "dens"
type = "proj" #or "slice"
axis = "xz"
zoom = False #boolean
zoom_cell = 226
annotate_vel_freq = 2
start_frame = 0
source_directory = "/Users/rajikak/Output/omega_t_ff_0.4/CircumbinaryOutFlow_0.25_lref_10/"
directory = source_directory + "WIND_" + type + "_" + axis + "_000000"
if zoom == False:
    zoom_cell=0
    annotate_vel_freq = 8
year = 31557600
au = 1.496e13
m_times = []
its = []

#load times:
with open('movie_times.txt', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        m_times.append(float(row[0]))
print "loaded movie frame times"

#load initial file:
f = h5py.File(directory, 'r')
xmin = f['minmax_xyz'][0][0]/au
xmax = f['minmax_xyz'][0][1]/au
if axis == "xy":
    ymin = f['minmax_xyz'][1][0]/au
    ymax = f['minmax_xyz'][1][1]/au
else:
    ymin = f['minmax_xyz'][2][0]/au
    ymax = f['minmax_xyz'][2][1]/au
dim = np.shape(f[field + "_" + type + "_" + axis])[0]
cl = (xmax-xmin)/(dim)
x = np.arange(xmin+(zoom_cell*cl), xmax-(zoom_cell)*cl, cl)
y = np.arange(ymin+(zoom_cell*cl), ymax-(zoom_cell)*cl, cl)
x_vel = np.arange(xmin+(((annotate_vel_freq/2.)+zoom_cell)*cl), xmax-(zoom_cell*cl), cl*annotate_vel_freq)
y_vel = np.arange(ymin+(((annotate_vel_freq/2.)+zoom_cell)*cl), ymax-(zoom_cell*cl), cl*annotate_vel_freq)
X, Y = np.meshgrid(x, y)
X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
print "created meshs"

#find usable plots:
it = 0
mit = 0
diff_prev = np.inf

while mit < len(m_times)-1:
    mtime = m_times[mit]
    while it < max_file:
        directory = source_directory + "WIND_" + type + "_" + axis + "_" + ("%06d" % it)
        f = h5py.File(directory, 'r')
        diff = abs(mtime - float(str(f['time'][0])))
        if diff == 0:
            its.append(it)
            print "found", mtime
            mit = mit + 1
            if mit == len(m_times):
                break
            else:
                mtime = m_times[mit]
        elif np.isinf(diff_prev) and diff != 0:
            diff_prev = diff
        elif diff > diff_prev:
            its.append(it-1)
            print "found", mtime
            diff_prev = np.inf
            it = it - 1
            mit = mit + 1
            if mit == len(m_times):
                break
            else:
                mtime = m_times[mit]
        elif diff < diff_prev:
            diff_prev = diff
        if it == max_file-1 and mtime == m_times[-1]:
            its.append(it)
            print "found", mtime
        it = it + 1
print "found usable movie plots"

#find extremes for colour bar
max = []
min = []
for it in its:
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + ("%06d" % it)
    f = h5py.File(directory, 'r')
    if zoom:
        for x in range(zoom_cell, dim-zoom_cell):
            den_field = f[field + '_' + type + '_' + axis][x]
            if shape(den_field) == (dim,1):
                den_field = den_field.transpose()[0]
            if type == "proj":
                if axis == "xz":
                    den_field = den_field/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                else:
                    den_field = den_field/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
            max_val = np.max(den_field[zoom_cell:dim-zoom_cell])
            min_val = np.min(den_field[zoom_cell:dim-zoom_cell])
    else:
        max_val = np.max(f[field + '_' + type + '_' + axis]/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0]))
        min_val = np.min(f[field + '_' + type + '_' + axis]/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0]))
    max.append(max_val)
    min.append(min_val)
max_str = '%.1e' % np.max(max)
min_str = '%.1e' % np.min(min)
max_round = str(0.5*floor(2.0 * float(max_str[0:3])))
min_round = str(0.5*ceil(2.0 * float(min_str[0:3])))
cbar_max = float(max_round+max_str[-4:])
cbar_min = float(min_round+min_str[-4:])
print "found colour bar extremes"

#plot figures
frame_val = start_frame
its = its[start_frame:]
for it in its:
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + ("%06d" % it)
    f = h5py.File(directory, 'r')
    velx = []
    vely = []
    image = []
    magx = []
    magy = []
    for x in range(zoom_cell, dim-zoom_cell):
        image_val = f[field + "_" + type + "_" + axis][x]
        magx_val = f["mag"+axis[0]+ "_" + type + "_"+axis][x]
        magy_val = f["mag"+axis[1]+ "_" + type + "_"+axis][x]
        if shape(image_val) == (dim,1):
            image_val = image_val.transpose()
            magx_val = magx_val.transpose()
            magy_val = magy_val.transpose()
        image_val = image_val[0][zoom_cell: dim-zoom_cell]/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
        magx_val = magx_val[0][zoom_cell: dim-zoom_cell]
        magy_val = magy_val[0][zoom_cell: dim-zoom_cell]
        image.append(image_val)
        magx.append(magx_val)
        magy.append(magy_val)
    image = np.array(image)
    magx = np.array(magx)
    magy = np.array(magy)
    for x in range(zoom_cell, dim-zoom_cell, annotate_vel_freq):
        velx_temp = f["vel" + axis[0] + "_" + type + "_" + axis][x]
        vely_temp = f["vel" + axis[1] + "_" + type + "_" + axis][x]
        if shape(velx_temp) == (dim, 1):
            velx_temp = velx_temp.transpose()
            vely_temp = vely_temp.transpose()
        velx.append(velx_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq])
        vely.append(vely_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq])
    plt.clf()
    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max))
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    Q = quiver(X_vel, Y_vel, velx, vely)
    plt.streamplot(X, Y, magx, magy)
    if len(f.keys()) > 12:
        if axis == "xy":
            part_pos_1 = f["particlepositions"][0]/au
            part_pos_2 = f["particlepositions"][1]/au
        else:
            part_pos_1 = []
            part_pos_2 = []
            ypos_sort = np.sort(f["particlepositions"][1])
            for ypos in ypos_sort:
                for part_pos in range(len(f["particlepositions"][1])):
                    if f["particlepositions"][1][part_pos] == ypos:
                        part_pos_1.append(f["particlepositions"][0][part_pos]/au)
                        part_pos_2.append(f["particlepositions"][2][part_pos]/au)
        part_mass_max = np.max(f["particlemasses"])
        part_mass = f["particlemasses"]/part_mass_max
        ax.scatter(part_pos_1, part_pos_2, c=part_mass, cmap=mpl.cm.gray)
    ax.set_title('time='+str((f['time'][0])/year)+'years')
    ax.set_xlabel('x (AU)')
    if axis == "xy":
        ax.set_ylabel('y (AU)', labelpad=-20)
    else:
        ax.set_ylabel('z (AU)', labelpad=-20)
    if zoom:
        ax.set_xlim([ceil(X[0][0]/10)*10, floor(X[-1][-1]/10)*10])
        ax.set_ylim([ceil(Y[0][0]/10)*10, floor(Y[-1][-1]/10)*10])
    else:
        ax.set_xlim([-1336.8983957219252, 1331.6761363636206])
        ax.set_ylim([-1336.8983957219252, 1331.6761363636201])
    plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".png", bbox_inches='tight')
    frame_val = frame_val + 1
    print 'done frame_val =', frame_val