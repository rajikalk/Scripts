#Reduced the G100 512 full sink data to save just the  time and mass accretion of every sink for the first 150kyr

import numpy as np
import pickle
import sys

full_data = sys.argv[1]
outfile = sys.argv[2]

file_open = open(full_data, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(full_data, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()

reduced_data = {}
import pdb
pdb.set_trace()
reduced_data.update({'time': global_data['time'].T[0][:SFE_5_ind]})
reduced_data.update({'m': global_data['m'][:SFE_5_ind]})
reduced_data.update({'x': global_data['x'][:SFE_5_ind]})
reduced_data.update({'y': global_data['y'][:SFE_5_ind]})
reduced_data.update({'z': global_data['z'][:SFE_5_ind]})
reduced_data.update({'ux': global_data['ux'][:SFE_5_ind]})
reduced_data.update({'uy': global_data['uy'][:SFE_5_ind]})
reduced_data.update({'uz': global_data['uz'][:SFE_5_ind]})

file = open(outfile, 'wb')
pickle.dump((reduced_data), file)
file.close()
