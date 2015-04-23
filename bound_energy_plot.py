#Read csv for total energy

import csv
import matplotlib.pyplot as plt

time = []
bound = []
unbound = []
Total_E = []
Total_TE = []
Total_KPE = []
Total_GPE = []
Total_KE_g = []
Total_KE_p = []
Total_GPE_g = []
Total_GPE_pg = []
Total_GPE_pp = []
bm = []
ubm = []
TM = []
it = 0
Msun = 1.9891e+33

with open('energy_conservation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[0])
            time.append(time_val)
            bound_val = (row[1])
            bound.append(bound_val)
            unbound_val = (row[2])
            unbound.append(unbound_val)
            Total_E_val = abs(float(row[3]))
            Total_E.append(Total_E_val)
            Total_TE_val = float(row[4])
            Total_TE.append(Total_TE_val)
            Total_KPE_val = float(row[5])
            Total_KPE.append(Total_KPE_val)
            Total_GPE_val = abs(float(row[6]))
            Total_GPE.append(Total_GPE_val)
            Total_KE_g_val = float(row[7])
            Total_KE_g.append(Total_KE_g_val)
            Total_KE_p_val = float(row[8])
            Total_KE_p.append(Total_KE_p_val)
            Total_GPE_g_val = abs(float(row[9]))
            Total_GPE_g.append(Total_GPE_g_val)
            Total_GPE_pg_val = abs(float(row[10]))
            Total_GPE_pg.append(Total_GPE_pg_val)
            Total_GPE_pp_val = abs(float(row[11]))
            Total_GPE_pp.append(Total_GPE_pp_val)
            bound_mass = float(row[12])/Msun
            bm.append(bound_mass)
            unbound_mass = float(row[13])/Msun
            ubm.append(unbound_mass)
            tm = bound_mass + unbound_mass
            TM.append(tm)
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(18.5, 10.5))
plt.subplot(111)
plt.plot(time, bm, label = 'Bound mass')
plt.plot(time, ubm, label = 'Unbound mass')
plt.plot(time, TM, label = 'Total mass of gas in grid')
"""
plt.plot(time, Total_E, 'k-', label = 'Total Energy')
plt.plot(time, Total_TE, 'r-.', label = 'Total Thermal Energy (TE)')
plt.plot(time, Total_KPE, 'g-.', label = 'Total Kinetic Energy (KE)')
plt.plot(time, Total_GPE, 'b-.', label = 'Total aboslute Gravitational Potential Energy (GPE)')
plt.plot(time, Total_KE_g, 'gx', label = 'Total KE of gas')
plt.plot(time, Total_KE_p, 'go', label = 'Total KE of particles')
plt.plot(time, Total_GPE_g, 'bx', label = 'Total absolute GPE of gas')
plt.plot(time, Total_GPE_pg, 'b^', label = 'Total absolute GPE of particles on gas')
plt.plot(time, Total_GPE_pp, 'bo', label = 'Total absolute GPE of particles on each other')
"""
plt.legend()
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.semilogy()
#plt.ylim([1.e19, 1.e47])
#fig.tight_layout()
#plt.xlim([0, 1.5])
plt.title("Unbound ratio of total gas in the grid using total energy")
plt.xlabel("Time (years)")
plt.ylabel("Mass ($M_\odot$)")
plt.savefig('Total_energy_bound_unbound.png')

