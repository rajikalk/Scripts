import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import pickle

arc_file = glob.glob('reduced_arcs/arc_variability.pkl')
sky_file = glob.glob('reduced_data/sky_line_variability.pkl')

arc_file = open(arc_file[0], 'r')
MJD_arc, RV_arc, RV_arc_err = pickle.load(arc_file)

sky_file = open(sky_file[0], 'r')
MJD_sky, RV_sky, RV_sky_err = pickle.load(sky_file)

print "NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!"
plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(MJD_arc, RV_arc, yerr=RV_arc_err, label='arcs', color='m', fmt='o')
ax.errorbar(MJD_sky, RV_sky, yerr=RV_sky_err*0.1, label='sky', color='r', fmt='o')
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
#plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation_comp_err.png')

plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(MJD_arc, RV_arc, label='arcs', color='m', fmt='o')
ax.scatter(MJD_sky, RV_sky, label='sky', color='r', fmt='o')
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
#plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation_comp_no_err.png')



