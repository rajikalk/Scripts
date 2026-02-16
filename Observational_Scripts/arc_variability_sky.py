import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import pickle

from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

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

dell_template = 0.1
wave_template=np.arange(90000)*dell_template + 3000

MJD = []
RV = []
RV_sig = []
RV_slits = []
RV_abs = []
RV_abs_sig = []
Exp_time = []
SNR = []
Sptype = []
sky_templates = []

if args.weight_direction != 'hor':
    slits = 38
else:
    slits = 12

templates = glob.glob('/Users/rajikak/tools/wifes_sky.fits')
abs_temps = glob.glob('/Users/rajikak/tools/absorption_spec_AF.fits')
abs_spec = pyfits.open(templates[0])

bad_intervals=([0,5500],[5700,5850],[6000,6100],[6700,6800])
abs_bad_ints =([0,5500],[5500,6860],[6910,7000])
#bad_intervals=([0,5500],[5700, 5850],[6000,6100],[6400, 6800])

good_file = np.ones(len(files))

for file in files:
    print("DOING FILE:", file)
    hdu = pyfits.open(file)
    MJD_obs = hdu[0].header['MJD-OBS']
    MJD.append(MJD_obs)
    exp_t = hdu[0].header['EXPTIME']
    obj = hdu[0].header['OBJECT']
    Exp_time.append(exp_t)
    flux_stamp,wave,sky = ps.read_and_find_star_p08(file, sky_rad=5, return_sky=True)#, manual_click=True)
    #Do this slit by slit:
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields("sptype")
    data = custom_simbad.query_object(obj)
    try:
        sptype = data['SP_TYPE'][-1]
    except:
        sptype = np.nan
    Sptype.append(sptype)
    
    '''
    for slit in range(12):
        if file == files[0]:
            RV_slits.append([])
        print "DOING SLIT:", slit
        spectrum,sig = ps.weighted_extract_spectrum(np.array([sky[slit]]))
        snr_val = np.max(spectrum)/np.median(spectrum)
        SNR.append(snr_val)
        if snr_val >= 3.0:
            rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, (bad_intervals), save_figures=False, save_dir='Images/', heliocentric_correction=0.0, fig_fn='Images/xcorr_slit_'+str(slit))
        else:
            rv = np.nan
        RV_slits[slit].append(rv)
    '''
    spectrum,sig = ps.weighted_extract_spectrum(sky)
    snr_val = np.max(spectrum)/np.median(spectrum)
    SNR.append(snr_val)
    if snr_val >= 5.0:
        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, (bad_intervals), save_figures=True, save_dir='Images/', heliocentric_correction=0.0, fig_fn='Images/xcorr_'+file.split('-')[-1].split('.')[0])
        '''
        spectrum_interp = np.interp(wave_template,wave[100:-100]*(1 - (0.0 - 0.0)/2.998e5),spectrum[100:-100])
        sky_templates.append(spectrum_interp)
        '''
        #ps.make_wifes_p08_template(file, 'templates/', rv=0.0, )
    else:
        print("SNR TOO LOW")
        rv = np.nan
        rv_sig = np.nan
    RV.append(rv)
    RV_sig.append(rv_sig)

    spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
    rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,abs_temps, (abs_bad_ints), save_figures=True, save_dir='Images/', heliocentric_correction=0.0, fig_fn='Images/xcorr_abs_'+file.split('-')[-1].split('.')[0])
    RV_abs.append(rv)
    RV_abs_sig.append(rv_sig)
    print('rv_abs =', rv)

'''
import pdb
pdb.set_trace()
sky_template = np.median(sky_templates, axis=0)
outfn = '/Users/rajikak/tools/wifes_sky.fits'
pyfits.writeto(outfn,sky_template,clobber=True)
'''

RV = np.array(RV)
RV_slits = np.array(RV_slits)

'''
RV_med = []
RV_std = []
for obs in range(len(RV_slits[0])):
    usable_ints = np.where(np.isnan(RV_slits[:,obs])==False)[0]
    RV_med.append(np.median(RV_slits[usable_ints,obs]))
    RV_std.append(np.std(RV_slits[usable_ints,obs])/len(usable_ints))
'''

#x_av = movingaverage(interval, r)
#plot(x_av, y)

print("NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!")

#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
'''
for slit in range(slits):
    ax.scatter(MJD, RV_slits[slit], marker='x') #, label=str(slit))
colormap = plt.cm.get_cmap('gist_rainbow') #nipy_spectral, Set1,Paired
colorst = [colormap(i) for i in np.linspace(0, 0.9,len(ax.collections))]
for t,j1 in enumerate(ax.collections):
    j1.set_color(colorst[t])
ax.errorbar(MJD, RV_med, yerr=RV_std,label='median over slits', fmt='o', color='b')
'''
#ax.scatter(MJD, RV_left_2, label='sky', color='r')
#ax.scatter(MJD, RV_abs, label='absorption lines', color='b')
plt.clf()
plt.errorbar(MJD, RV, yerr=RV_sig, fmt='o', color='r', label='sky')
plt.errorbar(MJD, RV_abs, yerr=RV_abs_sig, fmt='o', label='absorption', color='k')
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
#plt.xlim([MJD[0],MJD[-1]])
#plt.ylim([10,70])1
plt.legend(loc='best')
#plt.show()
plt.savefig('RV_variation_over_slits_using_skylines.png')

#WRITE OUT PICKLE FILE WITH RVs AND ERRORS
file = open('sky_line_variability.pkl', 'w+')
pickle.dump((MJD, RV), file)
file.close()