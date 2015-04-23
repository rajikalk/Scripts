import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

max_file = 2839
it = 0
field = "dens"
type = "proj" #or "slice"
axis = "xz"
zoom = False #boolean
source_directory = "/Users/rajikak/Output/CircumbinaryOutFlow_0.25_lref_8/"
directory = source_directory + "WIND_" + type + "_" + axis + "_000000"
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
axis_val=0
if axis == "xy":
    ymin = f['minmax_xyz'][1][0]/au
    ymax = f['minmax_xyz'][1][1]/au
else:
    ymin = f['minmax_xyz'][2][0]/au
    ymax = f['minmax_xyz'][2][1]/au
dim = np.shape(f[field + "_" + type + "_" + axis])[0]
cl = (xmax-xmin)/(dim)
if zoom:
    x = np.arange(xmin+(225*cl), xmax-(225)*cl, cl)
    y = np.arange(ymin+(225*cl), ymax-(225)*cl, cl)
    x_vel = np.arange(xmin+(229*cl), xmax-(225*cl), cl*2)
    y_vel = np.arange(ymin+(229*cl), ymax-(225*cl), cl*2)
x = np.arange(xmin, xmax, cl)
y = np.arange(ymin, ymax, cl)
X, Y = np.meshgrid(x, y)
x_vel = np.arange(xmin+(4*cl), xmax, cl*8)
y_vel = np.arange(ymin+(4*cl), ymax, cl*8)
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
            mit = mit + 1
            if mit == len(m_times):
                break
            else:
                mtime = m_times[mit]
        elif np.isinf(diff_prev) and diff != 0:
            diff_prev = diff
        elif diff > diff_prev:
            its.append(it-1)
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
        it = it + 1
print "found usable movie plots"

#find extremes for colour bar
max = []
min = []
for it in its:
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + ("%06d" % it)
    f = h5py.File(directory, 'r')
    max_val = np.max(f[field + '_' + type + '_' + axis])
    min_val = np.min(f[field + '_' + type + '_' + axis])
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
frame_val = 0
for it in its:
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + ("%06d" % it)
    f = h5py.File(directory, 'r')
    t = (f['time'][0])/year
    velx = []
    vely = []
    image = []
    magx = []
    magy = []
    for x in range(dim):
        if axis == "xy":
            image_temp = f[field + "_" + type + "_" + axis][x]
            image_val = image_temp.transpose()[0]
            magx_temp = f["mag"+axis[0]+ "_" + type + "_"+axis][x]
            magx_val = magx_temp.transpose()[0]
            magy_temp = f["mag"+axis[1]+ "_" + type + "_"+axis][x]
            magy_val = magy_temp.transpose()[0]
        else:
            image_val = f[field + "_" + type + "_" + axis][x][0]
            magx_val = f["mag"+axis[0]+ "_" + type + "_" +axis][x][0]
            magy_val = f["mag"+axis[1]+ "_" + type + "_" +axis][x][0]
        image.append(image_val)
        magx.append(magx_val)
        magy.append(magy_val)
    image = np.array(image)
    magx = np.array(magx)
    magy = np.array(magy)
    if zoom:
        for x in np.arange(225, dim-225, 2):
            velx.append(f["vel" + axis[0] + "_" + type + "_" + axis][x][0][225:-225:2])
            vely.append(f["vel" + axis[1] + "_" + type + "_" + axis][x][0][225:-225:2])
    else:
        for x in np.arange(0, dim, 8):
            velx.append(f["vel" + axis[0] + "_" + type + "_" + axis][x][0][0::8])
            vely.append(f["vel" + axis[1] + "_" + type + "_" + axis][x][0][0::8])
    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max))
    cbar = plt.colorbar(plot, pad=0.0)
    if type == "slice":
        cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    else:
        cbar.set_label('Projected Density ($g/cm^2$)', rotation=270, labelpad = 20)
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
    ax.set_title('time='+str(t)+'years')
    ax.set_xlabel('x (AU)')
    if axis == "xy":
        ax.set_ylabel('y (AU)')
    else:
        ax.set_ylabel('z (AU)')
    if zoom:
        ax.set_xlim([-150, 150])
        ax.set_ylim([-150, 150])
    else:
        ax.set_xlim([-1336.8983957219252, 1331.6761363636206])
        ax.set_ylim([-1336.8983957219252, 1331.6761363636201])
    plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".png")
    frame_val = frame_val + 1
    print 'done frame_val =', frame_val