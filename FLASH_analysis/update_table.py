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
        spin_val = Spin_dir.split('/')[-1].split('_')[-1]
        table_line = spin_val + '&' + Lref + '&'
        for Mach_val in Mach_vals:
            sink_evol_file = Spin_dir + '/' + Setup + '/' + Mach_val + '/' + Lref + '/sinks_evol.dat'
            with open(sink_evol_file, 'r') as f:
                print("reading ", sink_evol_file)
                reader = csv.reader(f, delimiter='  ')
                for row in reader:
                    import pdb
                    pdb.set_trace
        
