#!/usr/bin/env python
import h5py
import numpy as np
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import commands
from subprocess import call
import csv
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-z", "--zoom", help="Will movie be zoomed in?", default=False)
parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default = 4)
parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens")
parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
parser.add_argument("-if", "--init_frame", help="initial frame to start with", default = 0)
parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", default = 100)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

it = 0
field = args.field
if args.axis == "xz":
    axis = args.axis
    type = "proj"
else:
    axis = args.axis
    type = "slice"

path = sys.argv[1]
save_dir = sys.argv[2]
start_frame = int(args.init_frame)
source_directory = commands.getoutput('ls ' + path + 'WIND_' + type + '*').split('\n')
directory = source_directory[0]
lref = path.split('/')[-2].split('_')[-1]
ang_val = path.split('/')[-3].split('_')[-1]
movie_type = path.split('/')[-2].split('_')[0]
den_pert = path.split('/')[-2].split('_')[1]
year = 31557600
au = 1.496e13
Msun = 1.989e33
its = []
movie_dt = 2
r_acc = 0

#load initial file:
f = h5py.File(source_directory[-1], 'r')
max_time = float(str(f['time'][0]/year))
r_acc = f['r_accretion'][0]/au
if int(lref) < 10:
    line_rad = r_acc
else:
    box_len = (f['minmax_xyz'][0][1]*2)/au
    line_rad = 0.01*box_len
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
annotate_vel_freq = dim/32.
if args.zoom:
    zoom_times = args.zoom_times
    zoom_cell = ((1 - 1./zoom_times)/2.)*dim
else:
    zoom_times = 0
    zoom_cell = 0
annotate_vel_freq = int(((xmax-((zoom_cell)*cl)) - (xmin+((zoom_cell)*cl)))/32)
x = np.arange(xmin+((zoom_cell)*cl), xmax-((zoom_cell)*cl), cl)
y = np.arange(ymin+((zoom_cell)*cl), ymax-((zoom_cell)*cl), cl)
'''
x_vel = np.arange(xmin+(((annotate_vel_freq/2.)+zoom_cell)*cl), xmax-(zoom_cell*cl), cl*annotate_vel_freq)
y_vel = np.arange(ymin+(((annotate_vel_freq/2.)+zoom_cell)*cl), ymax-(zoom_cell*cl), cl*annotate_vel_freq)
'''
x_vel = np.arange(xmin+((zoom_cell)*cl), xmax-(zoom_cell*cl), cl*annotate_vel_freq)
y_vel = np.arange(ymin+((zoom_cell)*cl), ymax-(zoom_cell*cl), cl*annotate_vel_freq)

X, Y = np.meshgrid(x, y)
X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
print "created meshs"

#Find sink praticle creation time
sink_form_time = 0.0
for source in source_directory:
    f = h5py.File(source, 'r')
    if 'particlepositions' in f.keys():
        particles = True
        sink_form_time = float(str(f['time'][0]/year))
        print "sink formation time =", sink_form_time
        break
print "found sink particle creation time"

#generate movie times:
interval = 100./int(args.presink_frames)
presink = np.arange(1, 101, interval)
m_times = ((sink_form_time/np.log(presink[-1]))*np.log(presink))-sink_form_time
m_times = np.round(m_times/100)*100
m_times = np.array(m_times).tolist()
postsink = m_times[-1] + movie_dt
while postsink < (max_time-sink_form_time):
    m_times.append(postsink)
    postsink = postsink + movie_dt
print "created list of times"

#find usable plots:
it = 0
mit = 0
diff_prev = np.inf
mtime = m_times[mit]

while mit < len(m_times):
    while it < len(source_directory):
        f = h5py.File(source_directory[it], 'r')
        diff = abs(mtime-((f['time'][0]/year)-sink_form_time))
        if diff > diff_prev:
            its.append(it-1)
            it = it - 1
            print "found time", m_times[mit]
            mit = mit + 1
            if mit == len(m_times):
                break
            mtime = m_times[mit]
            diff_prev = np.inf
        else:
            diff_prev = diff
        it = it + 1
    if it == len(source_directory):
        mit = len(m_times)
'''
#find extremes for colour bar
max = []
min = []
v_min_temp = []
v_max_temp = []
for it in its:
    f = h5py.File(source_directory[it], 'r')
    if args.zoom:
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
    for x in range(zoom_cell, dim-zoom_cell, annotate_vel_freq):
        velx_temp = f["vel" + axis[0] + "_" + type + "_" + axis][x]
        vely_temp = f["vel" + axis[1] + "_" + type + "_" + axis][x]
        vel_temp = np.sqrt(np.square(velx_temp) + np.square(vely_temp))
        v_min_temp.append(np.min(vel_temp))
        v_max_temp.append(np.max(vel_temp))
    #max.append(max_val)
    #min.append(min_val)
    print "extremes found for", it
#max_str = '%.1e' % np.max(max)
#min_str = '%.1e' % np.min(min)
#max_round = str(0.5*floor(2.0 * float(max_str[0:3])))
#min_round = str(0.5*ceil(2.0 * float(min_str[0:3])))
v_min = np.min(v_min_temp)
v_max = np.max(v_max_temp)
v_scale = 0.05/v_max
#cbar_max = float(max_round+max_str[-4:])
#cbar_min = float(min_round+min_str[-4:])
#print "found colour bar extremes"
'''
cbar_max = 1.e-13
cbar_min = 1.e-17
print "found colour bar extremes"

#plot figures
frame_val = start_frame
its = its[start_frame:]
m_times = m_times[start_frame:]
for it in range(len(its)):
    f = h5py.File(source_directory[its[it]], 'r')
    velx = []
    vely = []
    image = []
    magx = []
    magy = []
    arrow_widths = []
    for x in range(int(zoom_cell), int(dim-zoom_cell)):
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
    for x in range(int(zoom_cell), int(dim-zoom_cell), int(annotate_vel_freq)):
        velx_temp = f["vel" + axis[0] + "_" + type + "_" + axis][x]
        vely_temp = f["vel" + axis[1] + "_" + type + "_" + axis][x]
        if shape(velx_temp) == (dim, 1):
            velx_temp = velx_temp.transpose()
            vely_temp = vely_temp.transpose()
        velx.append(velx_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq])
        vely.append(vely_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq])
        '''
        vel_temp = v_scale*(np.sqrt(np.square(velx_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq]) + np.square(vely_temp[0][zoom_cell:dim-zoom_cell:annotate_vel_freq])))
        for vel in vel_temp:
            arrow_widths.append(vel)
        '''

    # Create plots
    plt.clf()
    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Density ($g/cm^3$)', rotation=270, labelpad = 20)
    if frame_val > 0:
        plt.streamplot(X, Y, magx, magy, density=2, linewidth=0.25, minlength=0.5, arrowstyle='-')
    else:
        plt.streamplot(X, Y, magx, magy, density=2, linewidth=0.25, minlength=0.5)
    Q = quiver(X_vel, Y_vel, velx, vely, width=0.0003, minlength=0.0, scale=1.e7, headwidth=10, headlength=12)
    qk = quiverkey(Q, 0.9, 0.05, 500000, r'$5\rm{km}\rm{s}^{-1}$', labelpos='N', color='w', labelcolor='w',fontproperties={'weight': 'bold'}, labelsep=0.05)
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
        particle_text = ''
        for pos_it in range(len(part_pos_1)):
            ax.plot((part_pos_1[pos_it]-(line_rad), part_pos_1[pos_it]+(line_rad)), (part_pos_2[pos_it], part_pos_2[pos_it]), lw=1, color='k')
            ax.plot((part_pos_1[pos_it], part_pos_1[pos_it]), (part_pos_2[pos_it]-(line_rad), part_pos_2[pos_it]+(line_rad)), lw=1, color='k')
            ax.plot((part_pos_1[pos_it]-(r_acc), part_pos_1[pos_it]+(r_acc)), (part_pos_2[pos_it], part_pos_2[pos_it]), lw=0.5, c=part_color[pos_it])
            ax.plot((part_pos_1[pos_it], part_pos_1[pos_it]), (part_pos_2[pos_it]-(r_acc), part_pos_2[pos_it]+(r_acc)), lw=0.5, c=part_color[pos_it])
            circle = mpatches.Circle([part_pos_1[pos_it], part_pos_2[pos_it]], r_acc, fill=False, lw=1, edgecolor='k')
            ax.add_patch(circle)
            P_msun = str(np.round((part_mass[pos_it]*part_mass_max)/Msun, decimals=2))
            particle_text = particle_text+'$M_'+str(pos_it+1)+'$='+P_msun+' '
        plt.text((1./2.)*(xmax-((zoom_cell)*cl)), (xmin+((zoom_cell)*cl))-((1./5.)*(xmax-((zoom_cell)*cl))), r''+particle_text)
    time_val = 10.0*(np.floor(m_times[it]/10.0))
    ax.set_title('time = '+str(int(time_val))+' years')
    ax.set_xlabel('x (AU)')
    if axis == "xy":
        ax.set_ylabel('y (AU)', labelpad=-20)
    else:
        ax.set_ylabel('z (AU)', labelpad=-20)
    if args.zoom:
        ax.set_xlim([ceil(X[0][0]/10)*10, floor(X[-1][-1]/10)*10])
        ax.set_ylim([ceil(Y[0][0]/10)*10, floor(Y[-1][-1]/10)*10])
    else:
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([xmin, xmax])
    plt.savefig(save_dir + "movie_frame_" + ("%06d" % frame_val) + ".eps", format='eps', bbox_inches='tight')
    #plt.savefig("movie_frame_" + ("%06d" % frame_val) + ".pdf", format='pdf', bbox_inches='tight', dpi=600)
    print 'done frame_val =', frame_val
    frame_val = frame_val + 1
    plt.close()

#Convert files to jpeg.
file_list = commands.getoutput('ls '+ save_dir + '*.eps').split('\n')
file_counter = 1
for file in file_list:
    call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file, file[:-4]+'.jpg'])
    print "DONE", file_counter, "/", len(file_list)
    file_counter = file_counter + 1
    os.remove(file)

print "making movie"
movie_title = 'omega_t_ff_'+ ang_val + '_' + movie_type + '_' + den_pert + '_lref_' + lref + '.avi'
call(['python', '/home/100/rlk100/Scripts/make_movie.py', '-o', movie_title, save_dir+'*'])
print "movie made"
