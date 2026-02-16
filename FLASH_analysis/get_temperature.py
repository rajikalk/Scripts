#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np

files = sorted(glob.glob('*plt_cnt*'))
file = files[-1]
part_file = file[:10] + 'part' + file[-5:]

ds = yt.load(file, particle_filename=part_file)
dd = ds.all_data()
gamma = yt.YTArray(np.ones(np.shape(dd['dens'])), '')
gamma[np.where((dd['dens'].value >2.50e-16)&(dd['dens'].value<3.84e-13))[0]] = 1.1
gamma[np.where((dd['dens'].value >3.84e-13)&(dd['dens'].value<3.84e-8))[0]] = 1.4
gamma[np.where((dd['dens'].value >3.84e-8)&(dd['dens'].value<3.84e-3))[0]] = 1.1
gamma[np.where((dd['dens'].value >3.84e-3))[0]] = 5/3

K = (dd['pres']/yt.YTArray(np.power(dd['dens'].value,gamma.value), 'g/cm**3')).in_units('cm**2/s**2')

