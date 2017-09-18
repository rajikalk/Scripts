import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import pickle

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-w_dir", "--weight_direction", help="hor for along a slit and vert for across slits", default='hor')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

files = glob.glob('*.p08.fits')

MJD = []
RV = []
Exp_time = []

if args.weight_direction != 'hor':
    slits = 38
else:
    slits = 12

template = glob.glob('/Users/rajikak/tools/5*.fits')

bad_intervals=([0,5500],[5700, 5850],[6000,6100],[6700, 6800])
good_file = np.ones(len(files))

for file in files:
    print "DOING FILE:", file
    hdu = pyfits.open(file)
    MJD_obs = hdu[0].header['MJD-OBS']
    MJD.append(MJD_obs)
    Exp_time.append(hdu[0].header['EXPTIME'])
    flux_stamp,wave,sky = ps.read_and_find_star_p08(file, sky_rad=5, return_sky=True)#, manual_click=True)
    spectrum,sig = ps.weighted_extract_spectrum(sky)
    spectrum = np.sum(np.sum(sky,axis=0),axis=0)/len(np.where(sky[:,:,0]>0.0))
    plt.clf()
    plt.plot(wave, spectrum, label='spectrum')
    plt.xlim([5000, 7000])
    plt.savefig('Images/sky_spect' + file.split('-')[-1].split('.')[0] +'.png')
    rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,template, ([0,5400],[6870,6890]), save_figures=False, save_dir='Images/', heliocentric_correction=0.0, fig_fn='Images/xcorr_slit_'+file.split('-')[-1].split('.')[0])
    RV.append(rv)
    print 'rv =', rv
    
    '''
    for slit in range(slits):
        print "DOING SLIT:", slit
        #flux_stamp = sky[slit]
        #flux_stamp =  np.nan_to_num(flux_stamp)
        #flux_stamp = np.array([flux_stamp])
        #spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
        spectrum = (np.sum(sky,axis=1)[slit]/len(np.where(sky[slit,:,0]>0.0)[0]))
        plt.clf()
        #plt.plot(wave_template, spec_temp, label='template')
        plt.plot(wave, spectrum, label='spectrum')
        plt.xlim([5000, 7000])
        plt.savefig('Images/sky_spect_slit' + str(slit)+ '.png')
        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,template, ([0,5400],[6870,6890]), save_figures=False, save_dir='Images/', heliocentric_correction=0.0, fig_fn='Images/xcorr_slit_'+str(slit))
        RV[slit].append(rv)
        print 'rv =', rv
    '''

RV = np.array(RV)
'''
RV_med = []
RV_std = []
for obs in range(len(RV[0])):
    RV_med.append(np.median(RV[:,obs]))
    RV_std.append(np.std(RV[:,obs]))
'''

#x_av = movingaverage(interval, r)
#plot(x_av, y)

print "NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!"
plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(MJD, RV)
'''
for slit in range(slits):
    ax.scatter(MJD, RV[slit], label=str(slit))
'''
#plt.plot(MJD, RV_med, label='median RV')
#ax.errorbar(MJD, RV_med, yerr=RV_std,label='median RV')
#colormap = plt.cm.get_cmap('gist_rainbow') #nipy_spectral, Set1,Paired
#colorst = [colormap(i) for i in np.linspace(0, 0.9,len(ax.collections))]
#for t,j1 in enumerate(ax.collections):
#    j1.set_color(colorst[t])
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
plt.xlim([MJD[0],MJD[-1]])
plt.ylim([10,70])
#plt.legend(loc='best')
plt.savefig('Arc_variation_sky.png')

#WRITE OUT PICKLE FILE WITH RVs AND ERRORS
file = open('sky_line_variability.pkl', 'w+')
pickle.dump((MJD, RV), file)
file.close()