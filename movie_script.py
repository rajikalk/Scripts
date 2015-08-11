import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import csv
import commands


max_file = 2190
it = 0
field = "dens"
type = "proj" #or "slice"
axis = "xz"
zoom = False #boolean
zoom_cell = 123
annotate_vel_freq = 1
start_frame = 0
path = "/Users/rajikak/Output/omega_t_ff_0.2/CircumbinaryOutFlow_0.25_lref_10/"
source_directory = commands.getoutput('ls ' + path).split('\n')
directory = source_directory[0]
if zoom == False:
    zoom_cell=0
    annotate_vel_freq = 16
year = 31557600
au = 1.496e13
its = []
'''
#load times:
with open('movie_times.txt', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        m_times.append(float(row[0]))
print "loaded movie frame times"
'''
#load initial file:
f = h5py.File(path+source_directory[-1], 'r')
max_time = float(str(f['time'][0]/year))
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
x = np.arange(xmin+((zoom_cell)*cl), xmax-((zoom_cell)*cl), cl)
y = np.arange(ymin+((zoom_cell)*cl), ymax-((zoom_cell)*cl), cl)
x_vel = np.arange(xmin+(((annotate_vel_freq/2.)+zoom_cell)*cl), xmax-(zoom_cell*cl), cl*annotate_vel_freq)
y_vel = np.arange(ymin+(((annotate_vel_freq/2.)+zoom_cell)*cl), ymax-(zoom_cell*cl), cl*annotate_vel_freq)
X, Y = np.meshgrid(x, y)
X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
print "created meshs"

#Find sink praticle creation time
particles = False
sink_form_time = 0.0
while particles == False:
    for source in source_directory:
        directory = path + source
        f = h5py.File(directory, 'r')
        if 'particlepositions' in f.keys():
            particles = True
            sink_form_time = float(str(f['time'][0]/year))
            print "sink formation time =", sink_form_time
            break

#generate movie times:
presink = np.arange(1, 101, 4)
m_times = ((sink_form_time/np.log(presink[-1]))*np.log(presink))-sink_form_time
m_times = np.round(m_times/100)*100
m_times = np.array(m_times).tolist()
postsink = m_times[-1] + 5
while postsink < (max_time-sink_form_time):
    m_times.append(postsink)
    postsink = postsink + 5
print "created list of times"

#find usable plots:
it = 0
mit = 0
diff_prev = np.inf
'''
while mit < len(m_times)-1:
    mtime = m_times[mit] + max_time
    while it < max_file:
        directory = path + source_directory[it]
        f = h5py.File(directory, 'r')
        diff = abs(mtime - float(str(f['time'][0]/year)))
        #print diff
        if diff == 0:
            its.append(it)
            #print "found", mtime
            mit = mit + 1
            if mit == len(m_times):
                break
            else:
                mtime = m_times[mit] + max_time
        elif np.isinf(diff_prev) and diff != 0:
            diff_prev = diff
        elif diff > diff_prev:
            its.append(it-1)
            #print "found", mtime
            diff_prev = np.inf
            it = it - 1
            mit = mit + 1
            if mit == len(m_times):
                break
            else:
                mtime = m_times[mit] + max_time
        elif diff < diff_prev:
            diff_prev = diff
        if it == max_file-1 and mtime == m_times[-1]:
            its.append(it)
            #print "found", mtime
        it = it + 1
        print mtime
    mit = len(m_times)
print "found usable movie plots"

#find extremes for colour bar
max = []
min = []
for it in range(len(source_directory)):
    directory = path + source_directory[it]
    f = h5py.File(directory, 'r')
    if zoom:
        for x in range(zoom_cell, dim-zoom_cell):
            den_field = f[field + '_' + type + '_' + axis][x]
            if shape(den_field) == (1,dim):
                den_field = den_field.transpose()
            den_field = den_field[zoom_cell:dim-zoom_cell]
            max_val = np.max(den_field)
            min_val = np.min(den_field)
            if type == "proj":
                if axis == "xz":
                    max_val = max_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                    min_val = min_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                else:
                    max_val = max_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
                    min_val = min_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
    else:
        max_val = np.max(f[field + '_' + type + '_' + axis])
        min_val = np.min(f[field + '_' + type + '_' + axis])
        if type == "proj":
            if axis == "xz":
                max_val = max_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
                min_val = min_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
            else:
                max_val = max_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
                min_val = min_val/(f["minmax_xyz"][2][1]-f["minmax_xyz"][2][0])
    max.append(max_val)
    min.append(min_val)
    print "extremes found for", it
max_str = '%.1e' % np.max(max)
min_str = '%.1e' % np.min(min)
max_round = str(0.5*floor(2.0 * float(max_str[0:3])))
min_round = str(0.5*ceil(2.0 * float(min_str[0:3])))
cbar_max = float(max_round+max_str[-4:])
cbar_min = float(min_round+min_str[-4:])
print "found colour bar extremes"
'''
cbar_max = 1.e-12
cbar_min = 1.e-17
print "found colour bar extremes"

#plot figures
frame_val = start_frame
source_directory = source_directory[start_frame:]
for source in source_directory:
    directory = path + source
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
        image_val = image_val[0][zoom_cell: dim-zoom_cell]
        if type == "proj":
            image_val = image_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
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
    plot.set_rasterized(True)
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    res = plt.streamplot(X, Y, magx, magy, density=2, linewidth=1, minlength=0.9)
    lines = res.lines.get_paths()
    Q = quiver(X_vel, Y_vel, velx, vely)
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
        part_color = []
        for mass in part_mass:
            if mass/part_mass_max == 1:
                color = 'y'
            else:
                color = 'r'
            part_color.append(color)
        for pos_it in range(len(part_pos_1)):
            ax.plot((part_pos_1[pos_it]-(2.5*cl), part_pos_1[pos_it]+(2.5*cl)), (part_pos_2[pos_it], part_pos_2[pos_it]), lw=2, color='k')
            ax.plot((part_pos_1[pos_it], part_pos_1[pos_it]), (part_pos_2[pos_it]-(2.5*cl), part_pos_2[pos_it]+(2.5*cl)), lw=2, color='k')
            ax.plot((part_pos_1[pos_it]-(2.4*cl), part_pos_1[pos_it]+(2.4*cl)), (part_pos_2[pos_it], part_pos_2[pos_it]), lw=1, c=part_color[pos_it])
            ax.plot((part_pos_1[pos_it], part_pos_1[pos_it]), (part_pos_2[pos_it]-(2.4*cl), part_pos_2[pos_it]+(2.4*cl)), lw=1, c=part_color[pos_it])
            circle = mpatches.Circle([part_pos_1[pos_it], part_pos_2[pos_it]], 2.5*cl, fill=False, lw=3, edgecolor='k')
            ax.add_patch(circle)
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
    plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".eps", format='eps', bbox_inches='tight')
#plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".pdf", format='pdf', bbox_inches='tight', dpi=600)
    print 'done frame_val =', frame_val
    frame_val = frame_val + 1