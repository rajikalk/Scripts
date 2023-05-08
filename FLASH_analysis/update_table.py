import csv
import sys
import numpy as np
import pickle
import glob

Spins_dirs = sorted(glob.glob('/home/kuruwira/fast/Protostellar_spin/Flash_2023/*'))
Levels = ['Lref_9', 'Lref_10', 'Lref_11']
Mach_vals = ['Mach_0.1', 'Mach_0.2']
Setup = 'Binary'

for Lref in Levels:
    for Spin_dir in Spins_dirs:
        import pdb
        pdb.set_trace()
        
