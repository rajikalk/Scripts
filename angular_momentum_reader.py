#Read csv

import csv
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

time = []
L_tot = []
L_p = []
L_g = []
L_corr = []
L_err_l = []
L_err_h = []
L_pos = []
L_neg = []
x_p = []
y_p = []
z_p = []
x_g = []
y_g = []
z_g = []
x = []
y = []
z = []
com_sep = []
pp1x = []
pp1y = []
pp1z = []
pp2x = []
pp2y = []
pp2z = []
energy = []
E_plane = []
E_perp = []
cumul = []
cumul_val = 0.0
it = 0
lu = 10000000000000     #length units in cm
Rsun = 69600000000      #Radius of sun in cm

with open('angular_momentum_conservation_z.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if it != 0:
            time_val = float(row[0])
            time.append(time_val)
            L_tot_val = float(row[1])
            L_tot.append(L_tot_val)
            L_p_val = float(row[2])
            L_p.append(L_p_val)
            L_g_val = float(row[3])
            L_g.append(L_g_val)
            #L_corr_val = float(row[5])
            #L_corr.append(L_corr_val)
            '''
            err_l = float(row[4])
            L_err_l.append(err_l)
            err_h = float(row[5])
            L_err_h.append(err_h)
            L_pos_val = float(row[6])
            L_pos.append(L_pos_val)
            L_neg_val = float(row[7])*(-1.)
            L_neg.append(L_neg_val)
            x_p_val = float(row[1])
            x_p.append(x_p_val)
            y_p_val = float(row[2])
            y_p.append(y_p_val)
            z_p_val = float(row[3])
            z_p.append(z_p_val)
            x_g_val = float(row[4])
            x_g.append(x_g_val)
            y_g_val = float(row[5])
            y_g.append(y_g_val)
            z_g_val = float(row[6])
            z_g.append(z_g_val)
            x_val = float(row[7])
            x.append(x_val)
            y_val = float(row[8])
            y.append(y_val)
            z_val = float(row[9])
            z.append(z_val)
            com_val = ((float(row[10]))*lu)/Rsun
            com_sep.append(com_val)
            pp1x_val = float(row[4])
            pp1x.append(pp1x_val)
            pp1y_val = float(row[5])
            pp1y.append(pp1y_val)
            pp1z_val = float(row[6])
            pp1z.append(pp1z_val)
            pp2x_val = float(row[7])
            pp2x.append(pp2x_val)
            pp2y_val = float(row[8])
            pp2y.append(pp2y_val)
            pp2z_val = float(row[9])
            pp2z.append(pp2z_val)
            E_pl = float(row[15])
            E_plane.append(E_pl)
            E_pe = float(row[16])
            E_perp.append(E_pe)
            energy_val = float(row[17])
            energy.append(energy_val)
            cumul_val = float(row[18])
            cumul.append(cumul_val)
            '''
        else:
            it = 1
        

plt.clf()
fig = plt.figure(figsize=(8, 7))
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90, azim=90)
ax.plot(x_p, y_p, z_p, 'r-', label='CoM of particles')
ax.plot(x_g, y_g, z_g, 'b-', label='CoM of gas')
ax.plot(x, y, z, 'g-', label='Total CoM')
#ax.scatter(pp1x, pp1y, pp1z, 'g', label='particle 1')
#ax.scatter(pp2x, pp2y, pp2z, 'b', label='particle 2')
ax.legend()
ax.set_xlabel('X (gridunits)')
ax.set_ylabel('Y (gridunits)')
ax.set_zlabel('Z (gridunits)')
ax.set_xlim([0.48, 0.53])
ax.set_ylim([0.48, 0.53])
ax.set_zlim([0.48, 0.53])
plt.figure(figsize=(8, 7))
top = plt.subplot(211)
top.plot(time, energy, label='Total Energy Loss')
top.plot(time, E_plane, label = 'Orbital Plane')
top.plot(time, E_perp, label = 'Perpendicular')
top.legend(loc='upper right')
#top.set_title("Angular momentum loss per dump vs time")
#top.set_label("Time (years)")
top.set_label("Energy (erg)")
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
bot = plt.subplot(212)
bot.plot(time, cumul, label='cumulative energy loss')
#bot.set_title("Cumulative Angular Momentum loss")
bot.set_xlabel("Time (years)")
bot.set_ylabel("Energy (erg)")
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
plt.plot(time, L_tot, 'g-', label='Total $L_z$')
plt.plot(time, L_p, 'b:', label='Particle $L_z$')
plt.plot(time, L_g, 'r--', label='Gas $L_z$')
#plt.errorbar(time, L_tot, yerr=[L_err_l,L_err_h], label='Total $L_z$')
#plt.plot(time, L_corr, label='Total $L_z$ corrected for loss')
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
#plt.xlim([0.0, 0.5])
#plt.plot(time, L_pos, label='positive $L_z$')
#plt.plot(time, L_neg, label='negative $L_z$')
#plt.legend()
#plt.semilogy()
#plt.title("Angular momentum timeseries")
plt.xlabel("Time (years)")
plt.ylabel("Angular Momentum ($gcm^2s$)")
#plt.tight_layout()

plt.savefig('angular_momentum_timeseries_z.png')

