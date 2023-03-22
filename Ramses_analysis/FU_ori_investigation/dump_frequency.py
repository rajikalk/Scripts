import csv
import numpy as np
import matplotlib.pyplot as plt

file_data = "file_data.txt"
dump_time = []

with open(file_data, 'r') as data_file:
    reader = csv.reader(data_file, delimiter=' ')
    for row in reader:
        file_match = [s for s in row if 'output' in s]
        if len(file_match) > 0:
            time_ind = row.index(file_match[0]) - 1
            time_split = row[time_ind].split(':')
            time_dec = float(time_split[0] + str(int(time_split[1])/60)[1:])
            if len(dump_time) == 0:
                dump_time.append(time_dec)
            if dump_time[-1] != time_dec:
                dump_time.append(time_dec)
            
data_file.close()

dt = []
for time_it in range(1, len(dump_time)):
    if dump_time[time_it] < dump_time[time_it-1]:
        import pdb
        pdb.set_trace()
    else:
        dt_val = dump_time[time_it] - dump_time[time_it-1]
    dt.append(dt_val)

plt.clf()
plt.plot(dt)
plt.ylabel("time between dumps")
plt.savefig("dump_frequency.png")
