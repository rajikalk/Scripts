import csv
import sys
import numpy as np
import pickle

shell_out_file = sys.argv[1]

step_number = []
dt = []
step_chk = []
chk_number = []

with open(shell_out_file, 'r') as f:
    reader = csv.reader(f)
    line_counter = 0
    found_start = False
    line_counter = 0
    for row in reader:
        line_counter = line_counter + 1
        if len(row) > 0:
            if found_start == False and len(step_number) > 0:
                values = list(filter(None, row[0].split(' ')))
                try:
                    if int(values[0]) == next_step:
                        step_number.append(next_step)
                        dt.append(float(values[2]))
                        next_step = next_step + 1
                except:
                    pass
            if found_start:
                values = list(filter(None, row[0].split(' ')))
                curr_step = int(values[0])
                curr_dt = float(values[2])
                step_number.append(curr_step)
                dt.append(curr_dt)
                found_start = False
                next_step = curr_step + 1
            if row[0][:7] == '       ' and row[0][7] == 'n':
                found_start = True
            if '*** Wrote checkpoint file to SPIN_hdf5_chk_' in row[0]:
                curr_chk_number = int(row[0].split('*** Wrote checkpoint file to SPIN_hdf5_chk_')[-1].split('****')[0])
                if len(step_number) > 0:
                    step_chk.append(curr_step)
                    chk_number.append(curr_chk_number)

import matplotlib.pyplot as plt

plt.semilogy(step_number, dt)
for chk_it in range(len(step_chk)):
    plt.axvline(step_chk[chk_it])
    plt.text(step_chk[chk_it]+1, np.mean(dt), str(chk_number[chk_it]), rotation=90)
plt.xlabel('step number')
plt.ylabel('dt (s)')
plt.savefig(shell_out_file.split('/')[-1]+'.png')
