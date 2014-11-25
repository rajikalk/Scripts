#mother of all plots

import csv
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#define arrays:
time1 = []
time2 = []
time3 = []
time4 = []
sep1 = []
sep2 = []
sep3 = []
sep4 = []
mass1 = []
mass2 = []
mass3 = []
mass4 = []
den1 = []
den2 = []
den3 = []
den4 = []
pv1 = []
pv2 = []
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
L_tot = []
L_p = []
L_g = []
L_pos = []
L_neg = []
L_corr = []
it = 0
in_M = []
in_B = []
in_U = []
out_M = []
out_B = []
out_U = []
Total_M = []
Total_B = []
Total_U = []
PD_M = []
PD_B = []
PD_U = []

'''
with open('massfile.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[1])
            time1.append(time_val)
            Mi = float(row[11])
            in_M.append(Mi)
            #Mb = float(row[5])
            #in_B.append(Mb)
            #Mu = float(row[6])
            #in_U.append(Mu)
            Mo = float(row[13])
            out_M.append(Mo)
            #Mb = float(row[5])
            #out_B.append(Mb)
            #Mu = float(row[6])
            #out_U.append(Mu)
            TM_val = float(row[12])
            Total_M.append(TM_val)
            #TM_B = float(row[32])
            #Total_B.append(TM_B)
            #TM_U = float(row[33])
            #Total_U.append(TM_U)
            #dump_M = float(row[34])
            #PD_M.append(dump_M)
            #dump_B = float(row[35])
            #PD_B.append(dump_B)
            #dump_U = float(row[36])
            #PD_U.append(dump_U)
        else:
            it = 1

with open('separation_rotation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            t1 = float(row[0])
            time1.append(t1)
            s1 = float(row[1])
            sep1.append(s1)
            #d1 = float(row[2])
            #den1.append(d1)
            t2 = float(row[2])
            time2.append(t2)
            s2 = float(row[3])
            sep2.append(s2)
            #d2 = float(row[5])
            #den2.append(d2)
            #t3 = float(row[4])
            #time3.append(t3)
            #s3 = float(row[5])
            #sep3.append(s3)
            #d3 = float(row[8])
            #den3.append(d3)
            #t4 = float(row[6])
            #time4.append(t4)
            #s4 = float(row[7])
            #sep4.append(s4)
            #d4 = float(row[11])
            #den4.append(d4)
        else:
            it = 1
'''
it = 0
with open('energy_resolution.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = (row[0])
            time3.append(time_val)
            bound_val = (row[27])
            bound.append(bound_val)
            unbound_val = (row[28])
            unbound.append(unbound_val)
            tot_E_val = (row[29])
            tot_E.append(tot_E_val)
            tot_TE_val = (row[30])
            tot_TE.append(tot_TE_val)
            tot_KE_val = (row[31])
            tot_KE.append(tot_KE_val)
            tot_GPE_val = (row[32])
            tot_GPE.append(tot_GPE_val)
            tot_KE_g_val = (row[33])
            tot_KE_g.append(tot_KE_g_val)
            tot_KE_p_val = (row[34])
            tot_KE_p.append(tot_KE_p_val)
            tot_GPE_g_val = (row[35])
            tot_GPE_g.append(tot_GPE_g_val)
            tot_GPE_pg_val = (row[36])
            tot_GPE_pg.append(tot_GPE_pg_val)
            tot_GPE_pp_val = (row[37])
            tot_GPE_pp.append(tot_GPE_pp_val)
        else:
            it = 1
'''

it=0
with open('angular_momentum_conservation_z.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[0])
            time2.append(time_val)
            L_tot_val = float(row[1])
            L_tot.append(L_tot_val)
            L_p_val = float(row[2])
            L_p.append(L_p_val)
            L_g_val = float(row[3])
            L_g.append(L_g_val)
            L_corr_val = float(row[5])
            L_corr.append(L_corr_val)
            #L_neg_val = float(row[7])*(-1.)
            #L_neg.append(L_neg_val)
            print time_val, L_tot_val, L_p_val, L_g_val
	else:
            it = 1

f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
f.set_size_inches(8,7)
ax1.plot(time1, in_M, 'k', label='Total Mass in box')
ax1.plot(time1, in_B, 'r', label='Bound Mass')
ax1.plot(time1, in_U, 'g', label='Unbound Mass')
ax1.axvline(x=0.01, color='k')
ax1.axvline(x=0.04, color='k')
ax1.axvline(x=0.10, color='k')
ax1.axvline(x=0.06, color='b', ls='--')
ax1.axvline(x=0.02, color='r', ls=':')
ax1.axvline(x=0.03, color='r', ls=':')
ax1.axvline(x=0.05, color='r', ls=':')
ax1.axvline(x=0.07, color='r', ls=':')
ax1.axvline(x=0.08, color='r', ls=':')
ax1.axvline(x=0.09, color='r', ls=':')
ax1.annotate('(i)', xy=(0.47, 0.003), xycoords='data')
ax1.set_ylabel('$\Delta$Mass ($M_\odot$)')
ax2.plot(time1, out_M, 'k', label='Total Mass out box')
ax2.plot(time1, out_B, 'r', label='Bound Mass')
ax2.plot(time1, out_U, 'g', label='Unbound Mass')
ax2.axvline(x=0.01, color='k')
ax2.axvline(x=0.04, color='k')
ax2.axvline(x=0.10, color='k')
ax2.axvline(x=0.06, color='b', ls='--')
ax2.axvline(x=0.02, color='r', ls=':')
ax2.axvline(x=0.03, color='r', ls=':')
ax2.axvline(x=0.05, color='r', ls=':')
ax2.axvline(x=0.07, color='r', ls=':')
ax2.axvline(x=0.08, color='r', ls=':')
ax2.axvline(x=0.09, color='r', ls=':')
ax2.annotate('(ii)', xy=(0.47, -0.003), xycoords='data')
#ax2.set_ylim([0.0, 0.6])
ax2.set_ylabel('$\Delta$Mass ($M_\odot$)')
ax3.plot(time1, Total_M, 'k', label='Total Mass')
ax3.plot(time1, Total_B, 'r', label='Bound Mass')
ax3.plot(time1, Total_U, 'g', label='Unbound Mass')
ax3.axvline(x=0.01, color='k')
ax3.axvline(x=0.04, color='k')
ax3.axvline(x=0.10, color='k')
ax3.axvline(x=0.06, color='b', ls='--')
ax3.axvline(x=0.02, color='r', ls=':')
ax3.axvline(x=0.03, color='r', ls=':')
ax3.axvline(x=0.05, color='r', ls=':')
ax3.axvline(x=0.07, color='r', ls=':')
ax3.axvline(x=0.08, color='r', ls=':')
ax3.axvline(x=0.09, color='r', ls=':')
ax3.annotate('(iii)', xy=(0.47, 0.001), xycoords='data')
ax3.set_ylabel('$\Delta$Mass ($M_\odot$)')
ax4.plot(time1, PD_M, 'k', label='Total Mass')
ax4.plot(time1, PD_B, 'r', label='Bound Mass')
ax4.plot(time1, PD_U, 'g', label='Unbound Mass')
ax4.axvline(x=0.01, color='k')
ax4.axvline(x=0.04, color='k')
ax4.axvline(x=0.10, color='k')
ax4.axvline(x=0.06, color='b', ls='--')
ax4.axvline(x=0.02, color='r', ls=':')
ax4.axvline(x=0.03, color='r', ls=':')
ax4.axvline(x=0.05, color='r', ls=':')
ax4.axvline(x=0.07, color='r', ls=':')
ax4.axvline(x=0.08, color='r', ls=':')
ax4.axvline(x=0.09, color='r', ls=':')
ax4.annotate('(iv)', xy=(0.47, -0.0005), xycoords='data')
ax4.set_ylabel('$\Delta$Mass ($M_\odot$)')
ax4.set_xlabel('Time ($years$)')
ax4.set_xlim([0.0, 0.5])
f.subplots_adjust(hspace=0.12)

fig = plt.figure(1)
fig.clf()
fig.set_size_inches(8,8)
top = fig.add_subplot(311)
top.plot(time1, in_M, 'b', label='Total Mass in box')
#top.plot(time1, in_B, 'g', label='Bound Mass')
#top.plot(time1, in_U, 'r', label='Unbound Mass')
#top.legend()
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
top.set_ylabel('Time to Apastron ($years$)')
mid = fig.add_subplot(312)
mid.plot(time1, out_M, 'b', label='Total Mass out box')
#mid.plot(time1, out_B, 'g', label='Bound Mass')
#mid.plot(time1, out_U, 'r', label='Unbound Mass')
mid.set_ylabel('Orbital Period ($years$)')
#mid.set_xlabel('Time ($years$)')
#mid.set_ylim([0.0, 0.6])
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
bot = fig.add_subplot(313)
bot.plot(time1, Total_M, 'b', label='Total Mass')
#bot.plot(time1, Total_B, 'r', label='Bound Mass')
#bot.plot(time1, Total_U, 'g', label='Unbound Mass')
bot.set_ylabel('Eccentricity')
bot.set_xlabel('Time ($years$)')
#bot.set_ylim([0.0, 0.6])
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
#pd = fig.add_subplot(414)
#pd.plot(time1, PD_M, 'k', label='Total Mass')
#pd.plot(time1, PD_B, 'r', label='Bound Mass')
#pd.plot(time1, PD_U, 'g', label='Unbound Mass')
#pd.set_ylabel('|$\Delta$ mass| ($M_\odot$)')
#pd.set_xlabel('Time ($years$)')
#bot.set_ylim([0.0, 0.6])
#plt.axvline(x=0.01, color='k')
#plt.axvline(x=0.04, color='k')
#plt.axvline(x=0.10, color='k')
#plt.axvline(x=0.06, color='b', ls='--')
#plt.axvline(x=0.02, color='r', ls=':')
#plt.axvline(x=0.03, color='r', ls=':')
#plt.axvline(x=0.05, color='r', ls=':')
#plt.axvline(x=0.07, color='r', ls=':')
#plt.axvline(x=0.08, color='r', ls=':')
#plt.axvline(x=0.09, color='r', ls=':')

#fig.subplots_adjust(hspace = 0)
top = fig.add_subplot(411)
top.plot(time1, sep)
#top.set_ylim([0.0, 0.15])
#top.set_title('Separation')
top.set_ylabel('Separation ($R_\odot$)')
upper_middle = fig.add_subplot(412)
upper_middle.plot(time1, pv1)
upper_middle.plot(time1, pv2)
#middle.set_title('About orbital plane ('+str(lower_upper_bound)+'>z>'+str(upper_lower_bound)+')')
upper_middle.set_ylabel('Particle Velocity ($cms^{-1}$)')

lower_middle = fig.add_subplot(111)
lower_middle.set_ylim([0.e52, 7.e52])
lower_middle.set_xlim([0, 0.5])
lower_middle.plot(time2, L_tot, label='Total $L_z$')
lower_middle.plot(time2, L_p, label='Particle $L_z$')
lower_middle.plot(time2, L_g, label='Gas $L_z$')
#lower_middle.plot(time2, L_pos, label='positive $L_z$')
#lower_middle.plot(time2, L_neg, label='negative $L_z$')
#lower_middle.legend()
#middle.set_title('About orbital plane ('+str(lower_upper_bound)+'>z>'+str(upper_lower_bound)+')')
lower_middle.set_ylabel('Angular Momentum ($gcm^2s$)')
'''
fig = plt.figure()
fig.clf()
'''
plt.plot(time2, L_tot, 'k-', label = 'Total Angular Momentum')
plt.plot(time2, L_p, 'r:', label = 'Particle Momentum')
plt.plot(time2, L_g, 'g--', label = 'Gas Momentum')
plt.plot(time2, L_corr, 'k-.', label = 'Corrected Momentum')
'''
plt.plot(time3, tot_E, 'k-', label = 'Total Energy')
plt.plot(time3, tot_TE, 'b-', label = 'Total Thermal Energy')
plt.plot(time3, tot_KE, 'g-', label = 'Total Kinetic Energy')
plt.plot(time3, tot_GPE, 'r-', label = 'Total Potential Energy')
plt.plot(time3, tot_KE_g, 'g--', label = 'Total Kinetic Energy of gas')
plt.plot(time3, tot_KE_p, 'g:', label = 'Total Kinetic Energy of particles')
plt.plot(time3, tot_GPE_g, 'r--', label = 'Total Potential Energy of gas on gas')
plt.plot(time3, tot_GPE_pg, 'r-.', label = 'Total Potential Energy of particles on gas')
plt.plot(time3, tot_GPE_pp, 'r:', label = 'Total Potential Energy of particle on particle')

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

#bot.legend()
#bot.set_title('Below orbital plane (0.0>z>'+str(lower_upper_bound)+')')
plt.xlabel('Time ($years$)')
plt.ylabel('$L_z$ ($gcm^2s^{-1}$)')
#plt.ylim([-1.e47, 0.4e47])
'''
fig = plt.figure()
fig.clf()
#fig.set_size_inches(6,8)
plt.plot(time1, sep1, 'k-', label='nonrotating')
plt.plot(time2, sep2, 'b--' ,label='rotating')
#plt.plot(time3, sep3, 'r:', label='nonrot_highden')
#plt.plot(time4, sep4, label='175000K')

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
plt.axvline(x=0.06, color='b', ls='--')

#plt.legend(loc='upper right')
plt.xlabel('Time ($years$)')
plt.ylabel('Separation ($R_\odot$)')
#fig.tight_layout()
#plt.xlim([0.0, 0.5])
#plt.ylim([0.0, 5.e-6]
'''
filename = 'Angular_momentum_conservation.eps'
plt.savefig(filename)
