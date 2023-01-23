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
reduced_data.update({'time': global_data['time'].T[0]})
reduced_data.update({'m': global_data['m']})
reduced_data.update({'ux': global_data['ux']})

file = open(outfile, 'wb')
pickle.dump((reduced_data), file)
file.close()
