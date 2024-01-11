#!/usr/bin/env python
import sys
import yt
yt.enable_parallelism()
import glob
import my_flash_module as mym
import numpy as np

input_dir = sys.argv[1]
save_dir = sys.argv[2]
plot_time = yt.YTQuantity(4220, 'yr')

files = sorted(glob.glob(input_dir + '*plt_cnt*'))
fn = mym.find_files([plot_time.value], files)[0]
part_file = 'part'.join(fn.split('plt_cnt'))
ds = yt.load(fn, particle_filename=part_file)

dd = ds.all_data()

#Calculate separation vector
particle_mass = dd['particle_mass'][:2].in_units('msun')
com_x = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_posx'][:2].in_units('au'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com_y = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_posy'][:2].in_units('au'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com_z = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_posz'][:2].in_units('au'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com = yt.YTArray([com_x, com_y, com_z])

com_vel_x = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_velx'][:2].in_units('km/s'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com_vel_y = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_vely'][:2].in_units('km/s'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com_vel_z = np.sum(dd['particle_mass'][:2].in_units('msun') * dd['particle_velz'][:2].in_units('km/s'))/np.sum(dd['particle_mass'][:2].in_units('msun'))
com_vel = yt.YTArray([com_vel_x, com_vel_y, com_vel_z])


#Calculate orbital angular momentum vector mr x v
dx = dd['particle_posx'][:2].in_units('au') - com[0]
dy = dd['particle_posy'][:2].in_units('au') - com[1]
dz = dd['particle_posz'][:2].in_units('au') - com[2]
d_pos = yt.YTArray([dx, dy, dz]).T

dvx = dd['particle_velx'][:2].in_units('km/s') - com[0]
dvy = dd['particle_vely'][:2].in_units('km/s') - com[1]
dvz = dd['particle_velz'][:2].in_units('km/s') - com[2]
d_vel = yt.YTArray([dvx, dvy, dvz]).T

L_orb = dd['particle_mass'].value * np.cross(d_pos, d_vel).T

import pdb
pdb.set_trace()

