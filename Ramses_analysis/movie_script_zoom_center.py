#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import gc
import my_ramses_fields_short as myf
import csv

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
files = sorted(glob.glob(input_dir+"*/info*.txt"))
del input_dir
gc.collect()

#Get zoom_center from input
import pdb
pdb.set_trace()

