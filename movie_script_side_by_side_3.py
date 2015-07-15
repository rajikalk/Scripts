import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

max_file = [2363, 3067]
it = 0
field = "dens"
type = ["proj", "proj"] #or "slice"
axis = ["xz", "xz"]
zoom = [False, False] #boolean
zoom_cell = [0,0] #of 226
annotate_vel_freq = [8,8]
start_frame = 0
source_directory = ["/Users/rajikak/Output/SingleStarOutFlow_lref_9/", "/Users/rajikak/Output/CircumbinaryOutFlow_0.25_lref_9/"]
directory = [source_directory[0] + "WIND_" + type[0] + "_" + axis[0] + "_000000", source_directory[1] + "WIND_" + type[1] + "_" + axis[1] + "_000000"]
axes_titles = ["Single Star", "Binary Star"]
lref = 9
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
    while mit < len(m_times): #while the movie iterator (mit) is less than the length of the movie times.
        mtime = m_times[mit] #find current movie time
        while it < max_file[file]: #while it is less than all the movie files
            directory = source_directory[file] + "WIND_" + type[file] + "_" + axis[file] + "_" + ("%06d" % it)
            f = h5py.File(directory, 'r')
            diff = abs(mtime - float(str(f['time'][0]/year))) #find the different between the current movie time and file time.
            if np.isinf(diff_prev) and diff != 0: #you don't always get a perfect match, so find the time with the minimum difference. If previous difference is undefined. Define it.
                diff_prev = diff
            if diff == 0: #if its a perfect match, add it
                its[file].append(it)
                #print "found it =", it
                mit = mit + 1 #then go onto the next movie time to find.
                if mit == len(m_times): #unless we've finished the list.
                    break
                else:
                    mtime = m_times[mit] #we get the new movie time.
            elif diff > diff_prev: #Now we can start comparing times. Now if the current diff is larger than the previous diff, it means the previous file time is closer to the current movie time. So add the previous file.
                its[file].append(it-1)
                #print "found it =", it-1
                diff_prev = np.inf #then redefine the previous difference.
                it = it - 1 #now remember we added the previous file because its the closest to the current movie time. Now we want to know whether this file time is closest to the next movie time. So we jump back one.
                mit = mit + 1 #but we go onto the next movie time.
                if mit == len(m_times):
                    break
                else:
                    mtime = m_times[mit]
            elif diff < diff_prev: #if current difference is smaller than previous difference, it means we're still honing in on the closest time.
                diff_prev = diff #so update difference.
            if it == max_file[file]-1: #if we're reach the last file, add it.
                its[file].append(it)
                mit = len(m_times)
                it = max_file[file]
            elif mtime == m_times[-1]: #if we've reached the last movie time, add whatever file its up to.
                its[file].append(it)
                mit = len(m_times)
                it = max_file[file]
            it = it + 1
print "found usable movie plots"

#filter out plots to make both videos have same number of frames:
if len(its[0]) != len(its[1]):
    if len(its[0]) < len(its[1]):
        less_frames = 0
        more_frames = 1
    else:
        less_frames = 1
        more_frames = 0
    its[less_frames] = its[less_frames]
    its[more_frames] = its[more_frames][:len(its[less_frames])]

print "filtered files further"


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
    part_pos= []
    particles = []
    part_mass = []
    time = []
    for file in range(len(source_directory)):
        directory = source_directory[file] + "WIND_" + type[file] + "_" + axis[file] + "_" + ("%06d" % its[file][it])
        f = h5py.File(directory, 'r')
        time_val = f['time'][0]/year
        time.append(time_val)
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
            particles.append(True)
            ypos_sort = np.sort(f["particlepositions"][1])
            ypos_sort[:] = ypos_sort[::-1]
            for ypos in ypos_sort:
                for pos in range(len(f["particlepositions"][1])):
                    part_pos.append([])
                    if f["particlepositions"][1][pos] == ypos:
                        part_pos[file].append(f["particlepositions"][0]/au)
                        part_pos[file].append(f["particlepositions"][1]/au)
                        part_pos[file].append(f["particlepositions"][2]/au)
            part_mass_max = np.max(f["particlemasses"])
            part_mass.append(f["particlemasses"]/part_mass_max)
        else:
            particles.append(False)
            part_pos.append([])
            part_mass.append(0)
    plt.clf()
    fig, axes = plt.subplots(nrows=1,ncols=2)
    fig.subplots_adjust(wspace=0.0)
    fig.set_size_inches(15,6)
    axes[1].set_yticklabels(axes[1].get_yticklabels(), visible=False)
    for ax in range(len(axes)):
        plot = axes[ax].pcolormesh(X[ax], Y[ax], image[ax], cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max))
        axes[ax].streamplot(X[ax], Y[ax], magx[ax], magy[ax])
        Q = axes[ax].quiver(X_vel[ax], Y_vel[ax], velx[ax], vely[ax])
        if particles[ax]:
            if axis[ax] == "xz":
                axes[ax].scatter(part_pos[ax][0], part_pos[ax][2], c=part_mass[ax], cmap=mpl.cm.gray)
            else:
                axes[ax].scatter(part_pos[ax][0], part_pos[ax][1], c=part_mass[ax], cmap=mpl.cm.gray)
        axes[ax].set_xlabel('x (AU)')
        axes[ax].set_title(axes_titles[ax])
        if zoom[ax]:
            axes[ax].set_xlim([ceil(X[ax][0][0]/10)*10, floor(X[ax][-1][-1]/10)*10])
            axes[ax].set_ylim([ceil(Y[ax][0][0]/10)*10, floor(Y[ax][-1][-1]/10)*10])
        else:
            axes[ax].set_xlim([-width_au[ax], width_au[ax]])
            axes[ax].set_ylim([-width_au[ax], width_au[ax]])
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.8, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(plot, cax=cbar_ax)
    cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    fig.suptitle('lref=' + str(lref) + ',    time='+str(sum(time)/len(time))+'years')
    if axis[0] == "xy":
        axes[0].set_ylabel('y (AU)', labelpad=-20)
    else:
        axes[0].set_ylabel('z (AU)', labelpad=-20)
    plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".png", bbox_inches='tight')
    print 'done frame_val =', frame_val
    frame_val = frame_val + 1