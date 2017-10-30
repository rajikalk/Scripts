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
MJD_sky, RV_sky = pickle.load(sky_file)
RV_sky[np.where(RV_sky<0.0)[0]] = np.nan


MJD_full = np.concatenate((MJD_arc,MJD_sky))
RV_full = np.concatenate(((RV_arc+np.median(np.nan_to_num(RV_sky))),RV_sky))
RV_full[np.where(RV_full<0.0)[0]] = np.nan

RV_sort = RV_full[np.argsort(MJD_full)]
MJD_sort = MJD_full[np.argsort(MJD_full)]
RV_sort[np.where(np.nan_to_num(RV_sort)==0)[0]] = np.median(np.nan_to_num(RV_sky))

N = 3
cumsum, moving_aves = [0], [RV_sort[0]]

for i, x in enumerate(RV_sort, 1):
    cumsum.append(cumsum[i-1] + x)
    if i>=N:
        moving_ave = (cumsum[i] - cumsum[i-N])/N
        #can do stuff with moving_ave here
        moving_aves.append(moving_ave)

moving_aves.append(RV_sort[-1])

diff = RV_sort - moving_aves

bins = np.arange(np.floor(np.min(diff)), np.ceil(np.max(diff)))
hist, bins = np.histogram(diff, bins=bins)
width = (bins[1] - bins[0])
centers = bins[:-1] + width*0.5


print "NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!"
plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(MJD_arc, RV_arc + np.median(np.nan_to_num(RV_sky)), yerr=RV_arc_err, label='arcs', color='m', fmt='o')
ax.scatter(MJD_sky, RV_sky, label='sky', color='r')
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
#plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation_comp_err.png')

plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(MJD_arc, RV_arc + np.median(np.nan_to_num(RV_sky)), label='arcs', color='m')
ax.scatter(MJD_sky, RV_sky, label='sky', color='r')
ax.plot(MJD_sort, moving_aves)
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
#plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation_comp_no_err.png')

plt.clf()
plt.bar(centers, hist, align='center', width=width)
plt.xlabel('RV from moving average')
plt.savefig('histogram of RV spread.png')


