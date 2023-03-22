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

dt = np.array(dump_time)[1:] - np.array(dump_time)[:-1]
plt.clf()
plt.plot(dt)
plt.ylabel("time between dumps")
plt.savefig("dump_frequency.png")
