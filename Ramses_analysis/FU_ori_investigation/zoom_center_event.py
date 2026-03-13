#!/usr/bin/env python

import glob
import sys
import csv

sink_id = sys.argv[1]

fname = '/home/100/rlk100/gdata/RAMSES/Global/G100/512_Resolution/stars_red_512.pkl'
glob_data = open(fname, 'rb')
sink_data = pickle.load(glob_data)
glob_data.close()

file = sorted(glob.glob('data/*/info*txt'))[0]
with open(file, "r") as f:
    reader = csv.reader(f, delimiter=' ')
    for row in reader:
        if row[0] == 'time':
            curr_time = eval(row[-1])
            break

import pdb
pdb.set_trace()
 
