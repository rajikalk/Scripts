import h5py
import glob
import numpy as np

files = sorted(glob.glob('*plt_cnt*'))

field_names = np.array([['dens'],['pres'],['velx'],['vely'],['velz'],['gpot'],['magx'],['magy'],['magz'], ['alfl']])

for file in files:
    f = h5py.File(file, 'r+')
    del f['unknown names']
    f.create_dataset('unknown names', data=field_names, dtype="S10")
    f.close()
    print('updated file:', file)
