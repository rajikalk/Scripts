import csv
import numpy
import matplotlib.pyplot as plt

file_data = "file_data.txt"
dump_time = []

with open(file_data, 'r') as data_file:
    reader = csv.reader(data_file, delimiter=' ')
    for row in reader:
        if len(row) > 11:
            if 'output' in row[11]:
                time_split = row[10].split(':')
                time_dec = float(time_split[0] + str(int(time_split[1])/60)[1:])
                if len(dump_time) == 0:
                    dump_time.append(time_dec)
                if dump_time[-1] != time_dec:
                    dump_time.append(time_dec)
            
f.close()

dt = np.array(dump_time)[1:] - np.array(dump_time)[:-1]
plt.clf()
plt.plot(dt)
plt.ylabel("time between dumps")
plt.savefig("dump_frequency.png")
