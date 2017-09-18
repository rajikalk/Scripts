import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-w_dir", "--weight_direction", help="hor for along a slit and vert for across slits", default='hor')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

files = glob.glob('./*p08*')

MJD = []
RV = []
Amb_Temp = []
Baro_Pres = []
Rel_humid = []
CCD_Temp = []
Cam_Temp = []

if args.weight_direction != 'hor':
    slits = 38
else:
    slits = 12

#Make template of arc for each slitlet at the same pixel
print "MAKING TEMPLATES"
#read in file
temp_file = files[0]
hdu = pyfits.open(temp_file)
amb = hdu[0].header['TDK']
Amb_Temp.append(amb)
bar_pres = hdu[0].header['PMB']
Baro_Pres.append(bar_pres)
rel_hum = hdu[0].header['RH']
Rel_humid.append(rel_hum)
ccd = hdu[0].header['CCDTEMP']
CCD_Temp.append(ccd)
cam = hdu[0].header['CAMTEMPB']
Cam_Temp.append(cam)


#Get fluxes
flux = np.array([hdu[i].data for i in range(1,13)])
#Get wavelength scale and convert to logscale
wave_0 = hdu[1].header['CRVAL1'] + np.arange(flux.shape[2])*hdu[1].header['CDELT1']
dell_template = 0.1
wave_template=np.arange(90000)*dell_template + 3000
#get MJD to save
Obs_date = hdu[0].header['DATE-OBS'].split('T')[0]
MJD_obs = hdu[0].header['MJD-OBS']
MJD.append(MJD_obs)
for slit in range(slits):
    RV.append([0.0])
    if args.weight_direction != 'hor':
        flux_stamp = flux[:, slit, :]
    else:
        flux_stamp = flux[slit]
    flux_stamp =  np.nan_to_num(flux_stamp)
    flux_stamp = np.array([flux_stamp])
    spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
    dell_template = 0.1
    wave_template=np.arange(90000)*dell_template + 3000
    spectrum_interp = np.interp(wave_template,wave_0*(1 - (0.0 - 0.0)/2.998e5),spectrum)
    outfn = 'templates/slit_'+str(slit)+'.fits'
    pyfits.writeto(outfn,spectrum_interp,clobber=True)
    print('saved template:'+outfn)
print "MADE TEMPLATES FOR EACH SLIT"
print "NOW CALCULATING RV VARIATION FOR EACH SLIT"

bad_intervals= ([0,5400],[6700,7000]) #([0,5400],[6000,7000]) #([0,6100],[6700,7000])# ([0,5400],[6700,7000])

for file in files[1:]:
    print "DOING FILE:", file
    hdu = pyfits.open(file)
    flux = np.array([hdu[i].data for i in range(1,13)])
    #get wavelength scale and convert to log scale
    wave = hdu[1].header['CRVAL1'] + np.arange(flux.shape[2])*hdu[1].header['CDELT1']
    
    MJD_obs = hdu[0].header['MJD-OBS']
    MJD.append(MJD_obs)
    amb = hdu[0].header['TDK']
    Amb_Temp.append(amb)
    bar_pres = hdu[0].header['PMB']
    Baro_Pres.append(bar_pres)
    rel_hum = hdu[0].header['RH']
    Rel_humid.append(rel_hum)
    ccd = hdu[0].header['CCDTEMP']
    CCD_Temp.append(ccd)
    cam = hdu[0].header['CAMTEMPB']
    Cam_Temp.append(cam)

    for slit in range(slits):
        print "DOING SLIT:", slit
        if args.weight_direction != 'hor':
            flux_stamp = flux[:, slit, :]
        else:
            flux_stamp = flux[slit]
        flux_stamp =  np.nan_to_num(flux_stamp)
        flux_stamp = np.array([flux_stamp])
        spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
        template = ['templates/slit_'+str(slit)+'.fits']
        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,template, ([0,5400],[6870,6890]), save_figures=False, save_dir='Images/', heliocentric_correction=0.0)
        RV[slit].append(rv)
        print 'rv =', rv


RV = np.array(RV)
RV_med = []
RV_std = []
for obs in range(len(RV[0])):
    RV_med.append(np.median(RV[:,obs]))
    RV_std.append(np.std(RV[:,obs]))



print "NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!"
plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
for slit in range(slits):
    ax.scatter(MJD, RV[slit], label=str(slit))
#plt.plot(MJD, RV_med, label='median RV')
ax.errorbar(MJD, RV_med, yerr=RV_std,label='median RV')
colormap = plt.cm.get_cmap('gist_rainbow') #nipy_spectral, Set1,Paired
colorst = [colormap(i) for i in np.linspace(0, 0.9,len(ax.collections))]
for t,j1 in enumerate(ax.collections):
    j1.set_color(colorst[t])
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
plt.title(Obs_date)
#plt.ylim([-10.0, 10.0])
#plt.legend(loc='best')
plt.savefig('Arc_variation_Rel.png')
plt.clf()

RV = np.array(RV)
RV_med = []
for obs in range(len(RV[0])):
    RV_med.append(np.median(RV[:,obs]))

RV_med = np.array(RV_med)
Amb_Temp = np.array(Amb_Temp)
Rel_humid = np.array(Rel_humid)
CCD_Temp = np.array(CCD_Temp)
Cam_Temp = np.array(Cam_Temp)

plt.clf()
fig =plt.figure(figsize=(6, 10))
ax1 = plt.subplot(511)
plt.title(Obs_date)
plt.plot(MJD, RV_med)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylabel('median RV (km/s)')

ax2 = plt.subplot(512, sharex=ax1)
plt.plot(MJD, Amb_Temp)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.ylabel('ambient temp (K)')
'''
plt.subplot(5, 1, 3)
plt.plot(MJD, Baro_Pres)
plt.ylabel('Barometric pressure (hPa)')
'''
ax3 = plt.subplot(513, sharex=ax1)
plt.plot(MJD, Rel_humid)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.ylabel('relative humidity')

ax4 = plt.subplot(514, sharex=ax1)
plt.plot(MJD, CCD_Temp)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.ylabel('CCD temp (K)')

ax5 = plt.subplot(515, sharex=ax1)
plt.plot(MJD, Cam_Temp)
plt.ylabel('Cam temp (K)')

plt.xlabel('MJD')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('Instrument_values.png')

plt.clf()

MJD = np.array(MJD)
d_MJD = MJD[:-1] + (MJD[1:] - MJD[:-1])
d_RV_med = RV_med[1:] - RV_med[:-1]
d_Amb_Temp = Amb_Temp[1:] - Amb_Temp[:-1]
d_Rel_humid = Rel_humid[1:] - Rel_humid[:-1]
d_CCD_Temp = CCD_Temp[1:] - CCD_Temp[:-1]
d_Cam_Temp = Cam_Temp[1:] - Cam_Temp[:-1]

ax1 = plt.subplot(511)
plt.plot(d_MJD, d_RV_med)
plt.ylabel('$\Delta$ median RV (km/s)')

ax2 = plt.subplot(512, sharex=ax1)
plt.plot(d_MJD, d_Amb_Temp)
plt.ylabel('$\Delta$ ambient temp (K)')

ax3 = plt.subplot(513, sharex=ax1)
plt.plot(d_MJD, d_Rel_humid)
plt.ylabel('$\Delta$ relative humidity')

ax4 = plt.subplot(514, sharex=ax1)
plt.plot(d_MJD, d_CCD_Temp)
plt.ylabel('$\Delta$ CCD temp (K)')

ax5 = plt.subplot(515, sharex=ax1)
plt.plot(d_MJD, d_Cam_Temp)
plt.ylabel('Cam temp (K)')

plt.savefig('delta_values.png')

#WRITE OUT PICKLE FILE WITH RVs AND ERRORS
file = open('arc_variability.pkl', 'w+')
pickle.dump((MJD, RV_med, RV_std), file)
file.close()


