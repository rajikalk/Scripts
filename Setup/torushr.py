from yt.mods import *
import h5py
import numpy as np
import math
import os

# the directory for the data and hierarchy files to read in:

#datadir_particles = '/data/acomp/jstaff/enzo-fix/jc-version/256/0.6/old/DD0020/CE0020.cpu0000'
#rootdir = '/data/acomp/jstaff/enzo-fix/jc-version/256/0.6/old/DD0020/CE0020.hierarchy.hdf5'
#datadir_particles = '/data/acomp/jstaff/enzo-fix/qb/jc-version/256/0.6r3-2/DD0020/CE0020.cpu0000'
#rootdir = '/data/acomp/jstaff/enzo-fix/qb/jc-version/256/0.6r3-2/DD0020/CE0020.hierarchy.hdf5'
#datadir_particles = '/data/acomp/jstaff/enzo-fix/qb/jc-version/256/bigbox/DD0011/CE0011.cpu0000'
#datadir_particles = '/data/acomp/jstaff/enzo-fix/qb/agb_smaller/bigbox/DD0007/CE0007.cpu0000'
#rootdir = '/data/acomp/jstaff/enzo-fix/qb/agb_smaller/bigbox/DD0007/CE0007.hierarchy.hdf5'
#datadir_particles = '/data1/acomp/jstaff/enzo-fix/raijin/rajika/smallbox/DD0000/CE0000.cpu0000'
#rootdir = '/data1/acomp/jstaff/enzo-fix/raijin/rajika/smallbox/DD0000/CE0000.hierarchy.hdf5'
datadir_particles = '/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/highres/DD0000/CE0000.cpu0000'
rootdir = '/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/highres/DD0000/CE0000.hierarchy.hdf5'
#rootdir = '/data/acomp/jstaff/enzo-fix/qb/jc-version/256/bigbox/DD0011/CE0011.hierarchy.hdf5'
pf=load('/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/highres/DD0000/CE0000')
common_envelope = pf.h.all_data()
x=common_envelope['x']
y=common_envelope['y']
z=common_envelope['z']


# open hierarchy file and the subgrids:

hdf5_hierarchy = h5py.File(rootdir)
level0 = hdf5_hierarchy['Level0']
main_grid = level0['Grid00000001']
levellookuptable = hdf5_hierarchy['LevelLookupTable']

# read in number of particles variable and adjust it:
# (only thing that must be adjusted in this file)

numberofparticles=main_grid['NumberOfParticles']
#newnumberofparticles=numberofparticles.value+2

# open data file:

hdf5_particles_data = h5py.File(datadir_particles,'r+')
particles_main_grid = hdf5_particles_data['Grid00000001']
meta=hdf5_particles_data['Metadata']

# read in the fluid variables from the data file (leave these unchanged for
# now)

dark_matter_density = particles_main_grid['Dark_Matter_Density']
density = particles_main_grid['Density'].value
denhdf = particles_main_grid['Density']
grav_potential = particles_main_grid['Grav_Potential']
total_energy = particles_main_grid['TotalEnergy'].value
xvelocity = particles_main_grid['x-velocity'].value
yvelocity = particles_main_grid['y-velocity'].value
zvelocity = particles_main_grid['z-velocity']

# xvel = -0.72105 # Rob, 3Rstar separation
#xvel = -0.883229 # Rob, 2Rstar separtion
#xvel = -0.004829 # small AGB + 10 Mj, 1.5Rstar separtion, no gas velocity, set particle velocity further down
#boxsize=383 # 255

#xvelhdf[...]=xvelocity
# read in and change the particle variables:

particle_index=particles_main_grid['particle_index']
#new_particle_index = np.append(particle_index,1)

particle_position_x=particles_main_grid['particle_position_x'].value
particleposdtype = particle_position_x.dtype
#particle_position_x=np.append(particle_position_x,0.5) # 

particle_position_y=particles_main_grid['particle_position_y']
#particle_position_y=np.append(particle_position_y,0.57300) # <--- 1.1 stellar radii (1.1*1.326*10^13 cm) small AGB star, bigger box

particle_position_z=particles_main_grid['particle_position_z']
#particle_position_z=np.append(particle_position_z,0.5) # 

particle_velocity_x=particles_main_grid['particle_velocity_x']
#particle_velocity_x=np.append(particle_velocity_x,0.549647) # <--- 10 Mj at 1.1 R_agb (small agb) , biggerbox
#particle_velocity_x[0]=xvel     #


particle_velocity_y=particles_main_grid['particle_velocity_y']
#particle_velocity_y=np.append(particle_velocity_y,0.0) # 

particle_velocity_z=particles_main_grid['particle_velocity_z']
#particle_velocity_z=np.append(particle_velocity_z,0.0) # 

particle_type=particles_main_grid['particle_type']
#particle_type=np.append(particle_type,11)

#particle 2:
nparticles=2 # number of particles (to be used later when writing the files)
#new_particle_index = np.append(new_particle_index,2)
#particle_position_x=np.append(particle_position_x,0.5) # 
#particle_position_y=np.append(particle_position_y,0.6326) # <--- 2 stellar radii (2*1.326*10^13 cm) small AGB star, bigbox
#particle_position_z=np.append(particle_position_z,0.5) # 
#particle_velocity_x=np.append(particle_velocity_x,0.4077) # <--- 10 Mj at 2 R_agb (small agb) 2.3*10^6 cm/s, bigbox
#particle_velocity_y=np.append(particle_velocity_y,0.0) # 
#particle_velocity_z=np.append(particle_velocity_z,0.0) # 
#particle_type=np.append(particle_type,11)

# read the rest of the stuff from the hierarchy file:

gridend = main_grid['GridEndIndex']
gridstart = main_grid['GridStartIndex']
gridright = main_grid['GridRightEdge']
gridleft = main_grid['GridLeftEdge']
griddimension = main_grid['GridDimension']
gridglobalposition = main_grid['GridGlobalPosition']

# set the mass of the planet (mass density actually)
# this is just a fancy way of saying 1/255^3:
cellvolume = (gridright.value[0]-gridleft.value[0])*(gridright.value[1]-gridleft.value[1])*(gridright.value[2]-gridleft.value[2])/((gridend.value[0]-gridstart.value[0])*(gridend.value[1]-gridstart.value[1])*(gridend.value[2]-gridstart.value[2]))

#planetmass = 0.60 # 0.6 in solar masses
#planetmass = 0.009542 # 10 M_jupiter in solar masses
#planetdensity = planetmass/cellvolume
particle_mass=particles_main_grid['particle_mass'].value
#particle 1:
#new_particle_mass = np.append(particle_mass, planetdensity)
#particle 2:
#new_particle_mass = np.append(new_particle_mass, planetdensity)

particle_position_x=particle_position_x-1/sum(particle_mass*cellvolume)*sum(particle_mass*cellvolume*particle_position_x)+0.5

boxsize=384
Rtor = 0.2 # 1.e13 cm in a box that is 5.e13 cm
torr = 0.1 # 0.5e13 cm in a box that is 5.e13 cm

potparticle1=-1.661460e+03/(4*3.14)*particle_mass[0]*cellvolume/((x-particle_position_x[0])**2+(y-particle_position_y[0])**2+(z-particle_position_z[0])**2)**0.5
potparticle2=-1.661460e+03/(4*3.14)*particle_mass[1]*cellvolume/((x-particle_position_x[1])**2+(y-particle_position_y[1])**2+(z-particle_position_z[1])**2)**0.5
pot=potparticle1+potparticle2
# total_energy = total_energy/500.0 # about 200 K
#total_energy = total_energy/10.0 # about 17543 K
total_energy = total_energy/5.0 # about 35086 K
#total_energy = total_energy/1.0 # about 200,000K ?
for i in range(0,boxsize):
  for j in range(0,boxsize):
    for k in range(0,boxsize):
#      if ((x[i*256**2]-0.5)**2+(y[j*256]-0.5)**2+(z[k]-0.5)**2+Rtor**2-torr**2)**2-4*Rtor**2*((x[i*256**2]-0.5)**2+(y[j*256]-0.5)**2) < 0.0:
        density[k][j][i] = 0.5027399#/10.0 #corresponds to 1.0e-6 g/cc (/10)
        radius = ((x[i*boxsize**2]-0.5)**2+(y[j*boxsize]-0.5)**2)**0.5
        if (radius == 0.0): 
          radius = 0.001
# Kepler velocity is sqrt(-pot):
        yvelocity[k][j][i] = (x[i*boxsize**2]-0.5)/radius * (-pot[i*boxsize**2+j*boxsize+boxsize/2])**0.5*3.0/4.0
        xvelocity[k][j][i] = -(y[j*boxsize]-0.5)/radius * (-pot[i*boxsize**2+j*boxsize+boxsize/2])**0.5*3.0/4.0
#      else:
#        density[k][j][i] = 11.2488/10000.0
# create a new hierarchy file and put everything in there:

#f=h5py.File('/data/acomp/jstaff/enzo-fix/qb/jc-version/256/bigbox/10Mj1.1r/DD0011/CE0011.hierarchy.hdf5.new')
#f=h5py.File('/data/acomp/jstaff/enzo-fix/qb/agb_smaller/bigbox/DD0007/CE0007.hierarchy.hdf5.new')
#f=h5py.File('/data1/acomp/jstaff/enzo-fix/raijin/rajika/DD0000/CE0000.hierarchy.hdf5.new')
#f=h5py.File('/data/acomp/jstaff/enzo-fix/jc-version/256/0.6/old/DD0020/CE0020.hierarchy.hdf5.new')
#subgroup = f.create_group("Level0/Grid00000001")
#dset2=subgroup.create_dataset('GridDimension', griddimension.shape, griddimension.dtype, griddimension.value)
#dset3=subgroup.create_dataset('GridEndIndex', gridend.shape, gridend.dtype, gridend.value)
#dset4=subgroup.create_dataset('GridGlobalPosition',gridglobalposition.shape,gridglobalposition.dtype,gridglobalposition.value)
#dset5=subgroup.create_dataset('GridLeftEdge',gridleft.shape, gridleft.dtype, gridleft.value)
#dset6=subgroup.create_dataset('GridRightEdge',gridright.shape, gridright.dtype, gridright.value)
#dset7=subgroup.create_dataset('GridStartIndex',gridstart.shape, gridstart.dtype, gridstart.value)
#dset8=subgroup.create_dataset('NumberOfParticles', numberofparticles.shape, numberofparticles.dtype, newnumberofparticles) # this is the only new thing
#
#dset9=f.create_dataset('LevelLookupTable', levellookuptable.shape, levellookuptable. dtype, levellookuptable.value)
#f.close()

# create a new data file and put everything in there:

#f=h5py.File('/data/acomp/jstaff/enzo-fix/qb/jc-version/256/bigbox/10Mj1.1r/DD0011/CE0011.cpu0000.new')
#f=h5py.File('/data/acomp/jstaff/enzo-fix/qb/agb_smaller/bigbox/DD0007/CE0007.cpu0000.new')
f=h5py.File('/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/highres_hot2/CE0000.cpu0000.new')
#f=h5py.File('/data/acomp/jstaff/enzo-fix/jc-version/256/0.6/old/DD0020/CE0020.cpu0000.new')
subgroup = f.create_group("Grid00000001")
subgroup2 = f.create_group("Metadata")
subgroup2.attrs['LevelCycleCount'] = meta.attrs['LevelCycleCount']
dset2=subgroup.create_dataset('Dark_Matter_Density', dark_matter_density.shape, dark_matter_density.dtype, dark_matter_density.value)
dset3=subgroup.create_dataset('Density', density.shape, density.dtype, density)
dset4=subgroup.create_dataset('Grav_Potential', grav_potential.shape, grav_potential.dtype, grav_potential.value)
dset5=subgroup.create_dataset('TotalEnergy', total_energy.shape, total_energy.dtype,total_energy)
dset6=subgroup.create_dataset('particle_index', (nparticles,), particle_index.dtype, particle_index)
dset7=subgroup.create_dataset('particle_mass', (nparticles,), particle_mass.dtype, particle_mass )
dset8=subgroup.create_dataset('particle_position_x', (nparticles,), "float64", particle_position_x)
dset9=subgroup.create_dataset('particle_position_y', (nparticles,), "float64", particle_position_y)
dset10=subgroup.create_dataset('particle_position_z', (nparticles,), "float64", particle_position_z)
dset11=subgroup.create_dataset('particle_type', (nparticles,), ">i8", particle_type)
dset12=subgroup.create_dataset('particle_velocity_x', (nparticles,), "float64", particle_velocity_x)
dset13=subgroup.create_dataset('particle_velocity_y', (nparticles,), "float64", particle_velocity_y)
dset14=subgroup.create_dataset('particle_velocity_z', (nparticles,), "float64", particle_velocity_z)
dset15=subgroup.create_dataset('x-velocity', xvelocity.shape, xvelocity.dtype, xvelocity)
dset16=subgroup.create_dataset('y-velocity', yvelocity.shape, yvelocity.dtype, yvelocity)
dset17=subgroup.create_dataset('z-velocity', zvelocity.shape, zvelocity.dtype, zvelocity)
f.close()

# done!
