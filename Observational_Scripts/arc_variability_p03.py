import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-w_dir", "--weight_direction", help="hor for along a slit and vert for across slits", default='hor')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

files = glob.glob('./reduced_data/*p03*')
files = files[10:]
file_p08 = glob.glob('./reduced_data/*p08*')[0]

slit_temp = []
MJD = []
RV = []

#Get template
hdu_p08 = pyfits.open(file_p08)

file_temp = files[0]
hdu_temp = pyfits.open(file_temp)
MJD_obs = hdu_temp[0].header['MJD-OBS']
MJD.append(MJD_obs)
#Make templates
for slit in range(11):
    RV.append([0.0])
    flux_stamp = np.array([hdu_temp[slit+1].data[:]])
    spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
    pix_arr = np.arange(0, len(spectrum), 0.1)
    spectrum_interp = np.interp(pix_arr,np.arange(0, len(spectrum)),spectrum)
    if slit == 0:
        spectrum_temp = spectrum_interp
    slit_temp.append(spectrum_interp)
    outfn = 'templates/slit_'+str(slit+1)+'.fits'
    pyfits.writeto(outfn,spectrum_interp,clobber=True)
    print('saved template:'+outfn)

#import pdb
#pdb.set_trace()

#Now go through the other files
bad_intervals=([0,5400],[6700,7000])

for file in files[1:]:
    print "DOING FILE:", file
    hdu = pyfits.open(file)
    #get wavelength scale and convert to log scale
    
    MJD_obs = hdu[0].header['MJD-OBS']
    MJD.append(MJD_obs)
    for slit in range(11):
        spectrum_temp = slit_temp[slit]
        flux_stamp = np.array([hdu[slit+1].data[:]])
        spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
        pix_arr = np.arange(0, len(spectrum), 0.1)
        spectrum_interp = np.interp(pix_arr,np.arange(0, len(spectrum)),spectrum)
        cor = np.correlate(spectrum_interp,spectrum_temp,'same')
        x_arr = pix_arr - pix_arr[-1]/2.
        pix_shift = x_arr[np.argmax(cor)]
        RV[slit].append(pix_shift)
        print 'rv =', pix_shift


print "NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!"
plt.clf()
#marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
fig = plt.figure()
ax = fig.add_subplot(111)
for slit in range(11):
    ax.scatter(MJD[:13]+MJD[14:], RV[slit][:13]+RV[slit][14:], label=str(slit))
colormap = plt.cm.get_cmap('gist_rainbow') #nipy_spectral, Set1,Paired
colorst = [colormap(i) for i in np.linspace(0, 0.9,len(ax.collections))]
for t,j1 in enumerate(ax.collections):
    j1.set_color(colorst[t])
plt.xlabel('MJD')
plt.ylabel('pixel')
#plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation.png')