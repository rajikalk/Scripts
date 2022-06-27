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

import pdb
pdb.set_trace()
