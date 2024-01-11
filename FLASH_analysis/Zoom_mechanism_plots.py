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

#Altnerative calculation
Primary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[0], dd['particle_posy'].in_units('au')[0], dd['particle_posz'].in_units('au')[0]])
Secondary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[1], dd['particle_posy'].in_units('au')[1], dd['particle_posz'].in_units('au')[1]])
d_pos = Secondary_pos - Primary_pos

Primary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[0], dd['particle_vely'].in_units('km/s')[0], dd['particle_velz'].in_units('km/s')[0]])
Secondary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[1], dd['particle_vely'].in_units('km/s')[1], dd['particle_velz'].in_units('km/s')[1]])
d_vel = Secondary_vel - Primary_vel
L_vec = np.cross(d_pos, d_vel).T
L_vec_norm = L_vec/np.sqrt(np.sum(L_vec**2))

slc = yt.OffAxisSlicePlot(ds, L_vec_norm, 'dens', center=Primary_pos, width=(200, 'au'), north_vector=[0, 1, 0])

import pdb
pdb.set_trace()

