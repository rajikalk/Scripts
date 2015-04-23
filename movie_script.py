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
it = 323
field = "dens"
type = "proj" #or "slice"
axis = "xz"
source_directory = "/Users/rajikak/Output/CircumbinaryOutFlow_0.25/"
directory = source_directory + "WIND_" + type + "_" + axis + "_000000"
year = 365.25*24*60*60
au = 1.496e13

#load initial file:
f = h5py.File(directory, 'r')
xmin = f['minmax_xyz'][0][0]/au
xmax = f['minmax_xyz'][0][1]/au
zmin = f['minmax_xyz'][2][0]/au
zmax = f['minmax_xyz'][2][1]/au
dim = np.shape(f[field + "_" + type + "_" + axis])[0]
cl = (xmax-xmin)/(dim)
x = np.arange(xmin, xmax, cl)
z = np.arange(zmin, zmax, cl)
X, Z = np.meshgrid(x, z)
x_vel = np.arange(xmin+(4*cl), xmax+(0*cl), cl*8)
z_vel = np.arange(zmin+(4*cl), zmax+(0*cl), cl*8)
X_vel, Z_vel = np.meshgrid(x_vel, z_vel)
max_val = []
min_val = []

while it < max_file:
    print "it =", it
    it_str = ("%06d" % it)
    directory = source_directory + "WIND_" + type + "_" + axis + "_" + it_str
    f = h5py.File(directory, 'r')
    '''
    max = np.max(f['dens_proj_xz'])
    min = np.min(f['dens_proj_xz'])
    max_val.append(max)
    min_val.append(min)
    '''
    time = (f['time'][0])/year
    image = []
    velx = []
    velz = []
    magx = []
    magz = []
    for x in range(dim):
        image.append(f[field + "_" + type + "_" + axis][x][0])
        magx.append(f["magx_proj_xz"][x][0])
        magz.append(f["magz_proj_xz"][x][0])
    image = np.array(image)
    magx = np.array(magx)
    magz = np.array(magz)
    mag = np.sqrt(magx*magx + magz*magz)
    lw = 2*(mag/mag.max())
    for x in np.arange(0, dim, 8):
        velx.append(f["velx_proj_xz"][x][0][0::8])
        velz.append(f["velz_proj_xz"][x][0][0::8])

#    matplotlib.rcParams['xtick.direction'] = 'out'
#    matplotlib.rcParams['ytick.direction'] = 'out'

    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Z, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-3, vmax=2.5e2))
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Projected Density ($g/cm^2$)', rotation=270, labelpad = 20)
    Q = quiver(X_vel, Z_vel, velx, velz)
    plt.streamplot(X, Z, magx, magz) #linewidth=lw)
    if len(f.keys()) > 12:
        part_pos = [f["particlepositions"][0]/au, f["particlepositions"][2]/au]
        ax.plot(part_pos[0], part_pos[1], marker='o', color='y')
        #ax.plot(0.0, 0.0, marker='o', color='y')
    ax.set_title('time='+str(time)+'years')
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('z (AU)')
    plt.savefig("frame_" + ("%06d" % it) + ".png")

    it = it +1

