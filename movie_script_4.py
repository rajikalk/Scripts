import h5py
import numpy as np
from numpy import ma
from pylab import *
import matplotlib
import matplotlib as cm
import matplotlib.pyplot as plt
import matplotlib.quiver
from matplotlib.colors import LogNorm

max_file = 906
it = 0
field = "dens"
type = "proj" #or "slice"
axis = "xy"
source_directory = "/Users/rajikak/Output/SingleStarOutFlow_lref_8/"
directory = source_directory + "WIND_" + type + "_" + axis + "_000000"
year = 365.25*24*60*60
au = 1.496e13

#load initial file:
f = h5py.File(directory, 'r')
xmin = f['minmax_xyz'][0][0]/au
xmax = f['minmax_xyz'][0][1]/au
axis_val=0
if axis == "xy":
    axis_val = 1
else:
    axis_val = 2
ymin = f['minmax_xyz'][axis_val][0]/au
ymax = f['minmax_xyz'][axis_val][1]/au
dim = np.shape(f[field + "_" + type + "_" + axis])[0]
cl = (xmax-xmin)/(dim)
x = np.arange(xmin, xmax, cl)
y = np.arange(ymin, ymax, cl)
X, Y = np.meshgrid(x, y)
x_vel = np.arange(xmin+(4*cl), xmax+(0*cl), cl*8)
y_vel = np.arange(ymin+(4*cl), ymax+(0*cl), cl*8)
X_vel, Y_vel = np.meshgrid(x_vel, y_vel)
max_val = []
min_val = []

while it < max_file:
    it_str = ("%06d" % it)
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + it_str
    f = h5py.File(directory, 'r')
    '''
    max = np.max(f['dens_proj_xy'])
    min = np.min(f['dens_proj_xy'])
    max_val.append(max)
    min_val.append(min)
    '''
    time = (f['time'][0])/year
    image = []
    velx = []
    vely = []
    magx = []
    magy = []
    for x in range(dim):
        image_val = []
        magx_val = []
        magy_val = []
        if axis == "xy":
            image_temp = f[field + "_" + type + "_" + axis][x]
            image_val = image_temp.transpose()[0]
            magx_temp = f["mag"+axis[0]+"_proj_"+axis][x]
            magx_val = magx_temp.transpose()[0]
            magy_temp = f["mag"+axis[1]+"_proj_"+axis][x]
            magy_val = magy_temp.transpose()[0]
        else:
            image_val = f[field + "_" + type + "_" + axis][x][0]
            magx_val = f["mag"+axis[0]+"_proj_"+axis][x][0]
            magy_val = f["mag"+axis[1]+"_proj_"+axis][x][0]
        image.append(image_val)
        magx.append(magx_val)
        magy.append(magy_val)
    image = np.array(image)
    magx = np.array(magx)
    magy = np.array(magy)
    mag = np.sqrt(magx*magx + magy*magy)
    lw = 2*(mag/mag.max())
    for x in np.arange(0, dim, 8):
        velx.append(f["vel"+axis[0]+"_proj_"+axis][x][0][0::8])
        vely.append(f["vel"+axis[1]+"_proj_"+axis][x][0][0::8])

    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-3, vmax=5.e1))
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Projected Density ($g/cm^2$)', rotation=270, labelpad = 20)
    Q = quiver(X_vel, Y_vel, velx, vely)
    plt.streamplot(X, Y, magx, magy)
    if len(f.keys()) > 12:
        part_pos = [f["particlepositions"][0]/au, f["particlepositions"][2]/au]
        ax.plot(part_pos[0], part_pos[1], marker='o', color='y')
    ax.set_title('time='+str(time)+'years')
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_xlim([-1336.8983957219252, 1331.6761363636206])
    ax.set_ylim([-1336.8983957219252, 1331.6761363636201])
    plt.savefig("frame_" + ("%06d" % it) + ".png")
    print 'done it =', it
    it = it +1

