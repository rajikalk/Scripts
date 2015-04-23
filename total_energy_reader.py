#Read csv for total energy

import csv
import matplotlib.pyplot as plt

time = []
bound = []
unbound = []
"""
Total_E = []
Total_TE = []
Total_KPE = []
Total_GPE = []
Total_KE_g = []
Total_KE_p = []
Total_GPE_g = []
Total_GPE_pg = []
Total_GPE_pp = []
"""
it = 0

with open('energy_conservation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[0])/1e46
            time.append(time_val)
            bound_val = (row[1])
            bound.append(bound_val)
            unbound_val = (row[2])
            unbound.append(unbound_val)
            """
            Total_E_val = float(row[3])/1e46
            Total_E.append(Total_E_val)
            Total_TE_val = float(row[4])/1e46
            Total_TE.append(Total_TE_val)
            Total_KPE_val = float(row[5])/1e46
            Total_KPE.append(Total_KPE_val)
            Total_GPE_val = float(row[6])/1e46
            Total_GPE.append(Total_GPE_val)
            Total_KE_g_val = float(row[7])/1e46
            Total_KE_g.append(Total_KE_g_val)
            Total_KE_p_val = float(row[8])/1e46
            Total_KE_p.append(Total_KE_p_val)
            Total_GPE_g_val = float(row[9])/1e46
            Total_GPE_g.append(Total_GPE_g_val)
            Total_GPE_pg_val = float(row[10])/1e46
            Total_GPE_pg.append(Total_GPE_pg_val)
            Total_GPE_pp_val = float(row[11])/1e46
            Total_GPE_pp.append(Total_GPE_pp_val)
            """
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(7, 5))
plt.subplot(111)
plt.plot(time, unbound, label = 'unbound')
#plt.plot(time, Total_E, label = 'Total Energy')
#plt.plot(time, Total_TE, label = 'Total Thermal Energy (KE)')
#plt.plot(time, Total_KPE, label = 'Total Kinetic Energy (TE)')
#plt.plot(time, Total_GPE, label = 'Total Gravitational Potential Energy (GPE)')
#plt.plot(time, Total_KE_g, label = 'KE of gas')
#plt.plot(time, Total_KE_p, label = 'KE of particles')
#plt.plot(time,Total_GPE_g, label = 'GPE of gas')
#plt.plot(time, Total_GPE_pg, label = 'GPE of particles on gas')
#plt.plot(time, Total_GPE_pp, label = 'GPE of particles on each other')
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.semilogy()
plt.title('unbound gas ratio over time')
plt.xlabel("Time (years)")
plt.ylabel("ratio (%)")
plt.savefig('timeseries/Energy_unbound_timeseries.png')

