#Read csv

import csv
import matplotlib.pyplot as plt

time = []
Ef = []
Ee = []
Vf = []
Ve = []
tm = []
it = 0
Msun = 1.9891e+33

with open('bound_multiple_method_mass.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = (row[0])
            time.append(time_val)
            Ef_val = float(row[1])
            Ef.append(Ef_val)
            Ee_val = float(row[2])
            Ee.append(Ee_val)
            Vf_val = float(row[3])
            Vf.append(Vf_val)
            Ve_val = float(row[4])
            Ve.append(Ve_val)
            total = float(row[5])
            tm.append(total)
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(18.5, 10.5))
top = plt.subplot(211)
top.plot(time, tm, label = 'Total Mass of gas in grid')
top.set_ylabel('Total Mass ($M_\odot$)')
bot = plt.subplot(212)
bot.plot(time, Ef, 'r-', label = 'Total Energy method, full grid')
bot.plot(time, Ee, 'r--', label = 'Total Energy method, edge grids')
bot.plot(time, Vf, 'b-', label = '$v/v_{escape}$, full grid')
bot.plot(time, Ve, 'b--', label = '$v/v_{escape}$, edge grids')
bot.legend()
#plt.semilogy()
#plt.title("Amount of unbound gas using different methods for $0.75v_{keplerian}$")
bot.set_xlabel("Time (years)")
bot.set_ylabel("Unbound mass (%)")
fig.tight_layout()
plt.savefig('Multiple_methods_mass_timeseries.png')

