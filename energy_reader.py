#Read csv

import csv
import matplotlib.pyplot as plt

time = []
bound = []
unbound = []
tot_E = []
tot_TE = []
tot_KE = []
tot_GPE = []
tot_KE_g = []
tot_KE_p = []
tot_GPE_g = []
tot_GPE_pg = []
tot_GPE_pp = []
TM = []
E_corr = []
it = 0

with open('energy_conservation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[0])
            time.append(time_val)
            '''
            bound_val = float(row[15])
            bound.append(bound_val)
            unbound_val = float(row[16])
            unbound.append(unbound_val)
            '''
            tot_E_val = float(row[3])
            tot_E.append(tot_E_val)
            tot_TE_val = float(row[4])
            tot_TE.append(tot_TE_val)
            tot_KE_val = float(row[5])
            tot_KE.append(tot_KE_val)
            tot_GPE_val = float(row[6])
            tot_GPE.append(tot_GPE_val)
            tot_KE_g_val = float(row[7])
            tot_KE_g.append(tot_KE_g_val)
            tot_KE_p_val = float(row[8])
            tot_KE_p.append(tot_KE_p_val)
            tot_GPE_g_val = float(row[9])
            tot_GPE_g.append(tot_GPE_g_val)
            tot_GPE_pg_val = float(row[10])
            tot_GPE_pg.append(tot_GPE_pg_val)
            tot_GPE_pp_val = float(row[11])
            tot_GPE_pp.append(tot_GPE_pp_val)
            '''
            TM_val = float(row[14])
            TM.append(TM_val)
            
            E_corr_val = float(row[17])
            E_corr.append(E_corr_val)
            '''
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(8, 7))
plt.subplot(111)
'''
plt.plot(time, TM, label = 'Total Mass')
plt.plot(time, bound, label = 'Bound Mass')
plt.plot(time, unbound, label = 'Unbound Mass')
plt.axvline(x=0.01, color='k')
plt.axvline(x=0.04, color='k')
plt.axvline(x=0.10, color='k')
plt.axvline(x=0.06, color='b', ls='--')
plt.axvline(x=0.02, color='r', ls=':')
plt.axvline(x=0.03, color='r', ls=':')
plt.axvline(x=0.05, color='r', ls=':')
plt.axvline(x=0.07, color='r', ls=':')
plt.axvline(x=0.08, color='r', ls=':')
plt.axvline(x=0.09, color='r', ls=':')
'''
plt.plot(time, tot_E, 'k-', label = 'Total Energy')
#plt.plot(time, E_corr, 'k--', label = 'Corrected Total for mass loss')
plt.plot(time, tot_TE, 'b-', label = 'Total Thermal Energy')
plt.plot(time, tot_KE, 'g-', label = 'Total Kinetic Energy')
plt.plot(time, tot_GPE, 'r-', label = 'Total Potential Energy')
plt.plot(time, tot_KE_g, 'g--', label = 'Total Kinetic Energy of gas')
plt.plot(time, tot_KE_p, 'g:', label = 'Total Kinetic Energy of particles')
plt.plot(time, tot_GPE_g, 'r--', label = 'Total Potential Energy of gas on gas')
plt.plot(time, tot_GPE_pg, 'r-.', label = 'Total Potential Energy of particles on gas')
plt.plot(time, tot_GPE_pp, 'r:', label = 'Total Potential Energy of particle on particle')
'''
plt.legend(loc='upper right')
#plt.semilogy()
#plt.title("unbound mass timeseries")
'''
plt.xlabel("Time (years)")
#plt.ylabel("Mass ($M_sun$)")
plt.ylabel("Energy ($ergs$)")

plt.savefig('energy_timeseries.png')

