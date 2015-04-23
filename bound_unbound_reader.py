#Read csv

import csv
import matplotlib.pyplot as plt

time = []
bound = []
unbound = []
total_mass = []
mass_loss = []
it = 0
Msun = 1.9891e+33       #grams/Msun

with open('energy_bound_edge.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = (row[0])
            time.append(time_val)
            bound_val = float(row[3])/Msun
            bound.append(bound_val)
            unbound_val = float(row[4])/Msun
            unbound.append(unbound_val)
            m_tot = bound_val + unbound_val
            total_mass.append(m_tot)
            mass_loss_val = sum(total_mass)
            mass_loss.append(mass_loss_val)
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(18.5, 10.5))
plt.subplot(111)
plt.plot(time, bound, label = 'Bound mass')
plt.plot(time, unbound, label = 'Unbound mass')
plt.plot(time, total_mass, label = 'Total mass')
plt.plot(time, mass_loss, label = 'Mass_loss')
plt.legend()
#plt.semilogy()
plt.xlabel("Time (years)")
plt.ylabel("Mass ($M_\odot$)")
plt.savefig('Energy_bound_edge_timeseries.png')

