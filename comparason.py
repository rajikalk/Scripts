#Read csv

import csv
import matplotlib.pyplot as plt

time_75k = []
sep_75k = []
time_5k = []
sep_5k = []
time = []
sep = []
it = 0

with open('separation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_75k_val = (row[0])
            time_75k.append(time_75k_val)
            sep_75k_val = (row[1])
            sep_75k.append(sep_75k_val)
            time_5k_val = (row[2])
            time_5k.append(time_5k_val)
            sep_5k_val = (row[3])
            sep_5k.append(sep_5k_val)
            time_val = (row[4])
            time.append(time_val)
            sep_val = (row[5])
            sep.append(sep_val)
        else:
            it = 1
        

plt.clf()
#fig = plt.figure(figsize=(18.5, 15))
plt.subplot(111)
plt.plot(time_75k, sep_75k, label = '$0.75v_{kep}$')
plt.plot(time_5k, sep_5k, label = '$0.5v_{kep}$')
plt.plot(time, sep, label = 'non-rotating')
plt.legend()
#plt.semilogy()
plt.xlim([0.00, 0.12])
plt.title("Separation of particles")
plt.xlabel("Time (years)")
plt.ylabel("Separation ($R_\odot$)")
plt.savefig('separation_timeseries.png')

