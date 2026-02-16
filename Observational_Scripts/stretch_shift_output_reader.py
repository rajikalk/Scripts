import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

rows = []
best_shifts = []
best_stretches = []

file = sys.argv[1]
plot_var = sys.argv[2]

with open('stretch_shift_output.txt', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        row_val = float(row[1])
        if row_val == 1:
            rows.append([])
            best_shifts.append([])
            best_stretches.append([])
        shift = float(row[3])
        if shift != 0.0:
            rows[-1].append(row_val)
            best_shifts[-1].append(shift)
            stretch = float(row[4][:-1])
            best_stretches[-1].append(stretch)

fig = plt.figure()
ax = plt.subplot(111)
if 'shift' in plot_var:
    for i in range(len(rows)):
        ax.plot(rows[i], best_shifts[i], label='slit '+str(i+1))
    ax.set_xlabel('y pixel')
    ax.set_xlim(xmin = 9)
    ax.set_ylim([4025, 4140])
    ax.set_ylabel('shift value')
    ax.legend(bbox_to_anchor=(1.1, 1.05))
else:
    for i in range(len(rows)):
        ax.plot(rows[i], best_stretches[i], label='slit '+str(i+1))
    ax.set_xlabel('y pixel')
    ax.set_xlim(xmin = 9)
    ax.set_ylim([0.999, 1.0025])
    ax.set_ylabel('stretch value')
    ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()


