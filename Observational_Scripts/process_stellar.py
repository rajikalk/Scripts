"""After analysis with WiFeS, this suite of routines can extract a star optimally
and calculate its radial velocity.

WARNING - this code is extremely rough in its initial commit

example lines of code...

Executing from the code directory, e.g. with Margaret's output directory:
rv_process_dir('PROCESSED_DATA_DIRECTORY', outdir =/priv/mulga1/mstream/wifes/wifes/tools')

fn = 'T2m3wr-20140617.144009-0167.p11.fits'
flux,sig,wave = read_and_find_star_p11(fn)

fn = 'T2m3wr-20140617.144009-0167.p08.fits'

*** 5 lines below to run todcor ***
%run process_stellar
#fn = 'C:/Users/Margaret/MPhil/rvs_p08files/T2m3wb-20141110.103724-0816.p08.fits'
#flux,wave = read_and_find_star_p08(fn)


files = glob.glob('C:/Users/Margaret/MPhil/rvs_p08files_TTHor_2016/*.fits')
base_path = 'C:/Users/Margaret/MPhil/rvs_p08files_TTHor_2016/'
for file in files:  
    file_num = int(file.split('.')[1].split('-')[1])
    flux,wave = read_and_find_star_p08(file)
    spectrum,sig = weighted_extract_spectrum(flux)
    helcor = pyfits.getheader(file)['RADVEL']
    wave_log, spect_int, model_spect = calc_rv_todcor(spectrum,wave,sig,\
        ['C:/Users/Margaret/MPhil/RV_templates/8500g40p00k2v50.txt',\
        'C:/Users/Margaret/MPhil/RV_templates/5000g35p00k2v50.txt'],\
        alpha=0.1,out_fn='rvs.txt',jd=file_num,return_fitted=True,\
        heliocentric_correction=helcor)

    plt.clf()
    plt.plot(wave_log, spect_int, label='Data')
    plt.plot(wave_log, model_spect, label='Model')
    plt.legend()
    plt.xlabel('Wavelength')

*** lines below test todcor ***
binspect,binwave,binsig=make_fake_binary(spectrum, wave, sig, ['RV_templates/9000g40p00k2v150.txt','RV_templates/5250g35p00k2v150.txt'],0.5,-200,+200)
calc_rv_todcor(binspect,binwave,binsig,['RV_templates/9000g40p00k2v150.txt','RV_templates/5250g35p00k2v150.txt'],alpha=0.5)

rv,rv_sig = calc_rv_template(spectrum,wave,sig,'template_conv', ([0,5400],[6870,6890]))
rv,rv_sig = calc_rv_template(spectrum,wave,sig,template_fns, ([0,5400],[6870,6890]))

** To test the flux normalization of todcor **
margaret_dir = '/Users/mireland/Dropbox/Margaret Masters/RV_templates_flux/'
margaret_dir = '/Users/mireland/python/pywifes/tools/margaret/'
template_fns = glob.glob(margaret_dir+'*txt')
margaret_dir = '/Users/mireland/data/wifes/rvs_p08files/'

#flux, wave = read_and_find_star_p08(margaret_dir + 'T2m3wb-20141110.113949-0831.p08.fits')
#flux, wave = read_and_find_star_p08(margaret_dir + 'T2m3wb-20141108.120006-0389.p08.fits')
flux, wave = read_and_find_star_p08(margaret_dir + 'T2m3wb-20141110.114748-0832.p08.fits')
spectrum, sig = weighted_extract_spectrum(flux)
dummy = calc_rv_todcor(spectrum, wave,sig, [template_fns[3], template_fns[1]], bad_intervals=[[0,3810], [5005, 5028]], alpha=0.25, plotit=True)
dummy = calc_rv_todcor(spectrum, wave,sig, [template_fns[3], template_fns[1]], bad_intervals=[[0,4200], [5067,5075],[5500,6000]], alpha=0.25, plotit=True)
dummy = calc_rv_todcor(spectrum, wave,sig, [template_fns[2], template_fns[1]], bad_intervals=[[0,4200], [5067,5075],[5500,6000]], alpha=0.25, plotit=True)

"""

from __future__ import print_function
try:
    import pyfits
except:
    import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize as op
import pdb
import glob
import pickle
from readcol import readcol
from mpl_toolkits.mplot3d import Axes3D
from astropy.modeling import models, fitting
import sep

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    
    # print 'x = %d, y = %d'%(
    #     ix, iy)
    
    # assign global variable to access outside of function
    global coords
    coords.append((ix, iy))
    
    # Disconnect after 2 clicks
    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close(1)
    return

coords = []

def read_and_find_star_p11(fn, manual_click=False, npix=7, subtract_sky=True,sky_rad=2):
    """Read in a cube and find the star.
    Return a postage stamp around the star and the coordinates
    within the stamp
    
    NB This didn't really work as the details of flux calibration doesn't easily 
    enable optimal extraction. 
    
    This function should probably be REMOVED.
    """
    a = pyfits.open(fn)
    #Assume Stellar mode.
    flux = a[0].data[:,:,13:]
    sig = a[1].data[:,:,13:]
    image = np.median(flux,axis=0)
    maxpx = np.unravel_index(np.argmax(image[1:-1,1:-1]),image[1:-1,1:-1].shape)
    maxpx = (maxpx[0]+1,maxpx[1]+1)
    plt.clf()
    plt.imshow(image,interpolation='nearest')
    plt.plot(maxpx[1],maxpx[0],'wx')
    if subtract_sky:
        xy = np.meshgrid(range(image.shape[1]),range(image.shape[0]))
        dist = np.sqrt((xy[0]-maxpx[1])**2.0 + (xy[1]-maxpx[0])**2.0)
        sky = np.where( (xy[0] > 0) & (xy[1] > 0) & 
            (xy[0] < image.shape[1]-1) & (xy[1] < image.shape[0]-1) &
            (dist > sky_rad) & (dist < image.shape[1]))
        for i in range(flux.shape[0]):
            flux[i,:,:] -= np.median(flux[i,sky[0],sky[1]])
    ymin = np.min([np.max([maxpx[0]-3,0]),image.shape[0]-npix])
    xmin = np.min([np.max([maxpx[1]-3,0]),image.shape[1]-npix])
    flux_stamp = flux[:,ymin:ymin+npix,xmin:xmin+npix]
    sig_stamp = sig[:,ymin:ymin+npix,xmin:xmin+npix]
    wave = a[0].header['CRVAL3'] + np.arange(flux.shape[0])*a[0].header['CDELT3']
    return flux_stamp,sig_stamp,wave
    
def read_and_find_star_p08(fn, manual_click=False, npix=7, subtract_sky=True,sky_rad=2,fig_fn='', return_sky=False, max_px=[], return_maxpx=False):
    """Read in a cube and find the star.
    Return a postage stamp around the star and the wavelength scale
    
    NB This didn't really work as the details of flux calibration doesn't easily 
    enable optimal extraction.
    
    Parameters
    ----------
    fn: string
        filename
    npix: int
        Number of pixels to extract
    """
    a = pyfits.open(fn)
    Obj_name = a[0].header['OBJECT']
    Obs_date = a[0].header['DATE-OBS'].split('T')[0]
    RA = a[0].header['RA']
    DEC = a[0].header['DEC']
    #Assume Stellar mode.
    flux = np.array([a[i].data for i in range(1,13)])
    var = np.array([a[i].data for i in range(26,38)])
    wave = a[1].header['CRVAL1'] + np.arange(flux.shape[2])*a[1].header['CDELT1']
    image = np.median(flux,axis=2)
    data = image[:,4:]
    data = data.copy(order='C')
    objects = sep.extract(data, 1.5)
    objects['x'] = objects['x'] + 4
    #!!! 1->7 is a HACK - because WiFeS seems to often fail on the edge pixels !!!
    plt.clf()
    global fig
    fig = plt.figure(1)
    plt.imshow(image,interpolation='nearest')
    if max_px != []:
        maxpx = (max_px[0],max_px[1])
    elif manual_click == True:
        global coords

        # Call click func
        global cid
        cid = fig.canvas.mpl_connect('button_press_event', onclick)

        plt.show()
        maxpx = (int(round(np.min([coords[0][1], coords[1][1]]))), int(round(np.min([coords[0][0], coords[1][0]]))))
        coords = []
        '''
        image_stamp = image[maxpx[0]-3:maxpx[0]+4][:, maxpx[1]-3:maxpx[1]+4]
        x_pos = np.arange(np.shape(image_stamp)[0])
        y_pos = np.arange(np.shape(image_stamp)[1])
        X, Y = np.meshgrid(y_pos, x_pos)
        compx = [((np.sum(Y*image_stamp))/np.sum(image_stamp))+maxpx[0]-3, ((np.sum(X*image_stamp))/np.sum(image_stamp))+maxpx[1]-3]
        '''
    else:
        maxpx = np.unravel_index(np.argmax(image[:,10:-10]),image[:,10:-10].shape)
        maxpx = (maxpx[0],maxpx[1]+10)
        '''
        x_pos = np.arange(np.shape(image[:,10:-10])[0])
        y_pos = np.arange(np.shape(image[:,10:-10])[1])
        X, Y = np.meshgrid(y_pos, x_pos)
        compx = [((np.sum(Y*image[:,10:-10]))/np.sum(image[:,10:-10])), ((np.sum(X*image[:,10:-10]))/np.sum(image[:,10:-10]))+10]
        '''
    plt.clf()
    plt.imshow(image,interpolation='nearest')
    plt.plot(maxpx[1],maxpx[0],'wx')
    for obj in range(len(objects['x'])):
        plt.plot(objects['x'][obj],objects['y'][obj],'w+')
    plt.colorbar(fraction=0.0155, pad=0.0)
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')
    #plt.show()
    plt.title(str(Obj_name) + '_' + str(Obs_date) + '_(' + str(RA) + ',' + str(DEC) + ')')
    
    xy = np.meshgrid(range(image.shape[1]),range(image.shape[0]))
    dist = np.sqrt((xy[0]-maxpx[1])**2.0 + (xy[1]-maxpx[0])**2.0)
    sky = (xy[0] >= 0) & (xy[1] >= 0) & (xy[0] < image.shape[1]) & (xy[1] < image.shape[0]) & (dist > sky_rad) & (dist < image.shape[1]) & (xy[0] > 7)
    for obj in range(len(objects['x'])):
        dist = np.sqrt((xy[0]-objects['x'][obj])**2.0 + (xy[1]-objects['y'][obj])**2.0)
        sky = sky * ((xy[0] >= 0) & (xy[1] >= 0) & (xy[0] < image.shape[1]) & (xy[1] < image.shape[0]) & (dist > sky_rad) & (dist < image.shape[1]) & (xy[0] > 7))
    sky_temp = np.where(sky)
    if subtract_sky:
        for i in range(flux.shape[2]):
            flux[:,:,i] -= np.median(flux[sky_temp[0],sky_temp[1],i])
    sky = np.dstack([sky]*np.shape(flux)[-1])
    sky = sky*flux
    ymin = np.min([np.max([maxpx[0]-npix//2,0]),image.shape[0]-npix])
    xmin = np.min([np.max([maxpx[1]-npix//2,0]),image.shape[1]-npix])
    flux_stamp = flux[ymin:ymin+npix,xmin:xmin+npix,:]
    var_stamp = var[ymin:ymin+npix,xmin:xmin+npix,:]
    #flux_stamp = flux[int(compx[0]):int(compx[0])+2,int(compx[1]):int(compx[1])+2,:]
    if len(fig_fn)>0:
        plt.savefig(fig_fn, bbox_inches='tight')
    if return_sky and return_maxpx:
        return flux_stamp,wave,var_stamp,sky,maxpx
    elif return_sky:
        return flux_stamp,wave,var_stamp,sky
    elif return_maxpx:
        return flux_stamp,wave,var_stamp,maxpx
    else:
        return flux_stamp,wave,var_stamp
    
def weighted_extract_spectrum(flux_stamp, variance, readout_var=11.0):
    """Optimally extract the spectrum based on a constant weighting
    
    Based on a p08 file axis ordering.
    
    Readout variance is roughly 11 in the p08 extracted spectra
    
    Parameters
    ----------
    flux_stamp: numpy array
        nx x ny x nwave IFU image as a function of wavelength
    
    readout_var: float (optional)
        Readout variance in extracted spectrum in DN.
    
    TODO: 
    1) Look for and remove bad pix/cosmic rays.
    2) Remove dodgy constant for readout_var.
    """
    #Find the median flux over all wavelengths, limiting to be >0
    if np.max(np.median(flux_stamp,axis=2)) > 0.0:
        flux_med = np.maximum(np.median(flux_stamp,axis=2),0)
        var_med = np.maximum(np.median(variance,axis=2),0)
    else:
        flux_stamp = flux_stamp - np.min(flux_stamp)
        flux_med = np.maximum(np.median(flux_stamp,axis=2),0)
        var_med = np.maximum(np.median(variance,axis=2),0)
        #flux_med = np.median(flux_stamp,axis=2) - np.min(np.median(flux_stamp,axis=2))
    
    pixel_var = flux_med + readout_var
    weights = flux_med/pixel_var
    n_spaxels = np.prod(weights.shape)

    #Form a weighted average, then multiply by n_spaxels to get a sum
    spectrum = n_spaxels * np.array([np.sum(flux_stamp[:,:,i]*weights)/np.sum(weights) for i in range(flux_stamp.shape[2])]) 
    spectrum = np.array([np.sum(flux_stamp[:,:,i]*weights) for i in range(flux_stamp.shape[2])]) 
    
    #Old calculation of sigma.  Lets be a little more readable!
    sig = np.array([np.sqrt(np.sum((np.maximum(flux_stamp[:,:,i],0)+readout_var)*weights**2)) for i in range(flux_stamp.shape[2])])

    #The variance of each pixel is flux_stamp + readout_var, with flux_stamp being an estimate
    #of flux per pixel, which should not be less than zero.
    #var = [np.sum((np.maximum(flux_stamp[:,:,i],0)+readout_var)*weights**2)/np.sum(weights)**2 for i in range(flux_stamp.shape[2])]
    #sig = n_spaxels * np.sqrt(np.array(var))
    
    return spectrum,sig
    
def conv_ambre_spect(ambre_dir,ambre_conv_dir):
    """Take all the AMBRE spectra from a directory, convolve and re-sample
    by a factor of 10, then save to a new directory"""
    infns = glob.glob(ambre_dir + '/*fits')
    for infn in infns:
        data = pyfits.getdata(infn)
        data = np.convolve(data,np.ones(10)/10., 'same')
        conv_data = data[10*np.arange(90000,dtype='int')].astype('float32')
        ix_start = infn.rfind('/') + 1
        ix_end = infn.rfind('.')
        outfn = infn[ix_start:ix_end] + 'conv.fits'
        pyfits.writeto(ambre_conv_dir + '/' + outfn,conv_data, clobber=True)
    
def conv_tlusty_spect(tlusty_dir,tlusty_conv_dir):    
    """
    Take all Tlusty spectra from a directory, convole to 0.1A, 
    then save to a new directory
    Currently resampling onto a wavelength grid of 0.1A also, from 
    3000 to 12000A to match AMBRE spectra, note Tlusty only covers 3000A to 10000A
    also mostly matching filenames
    """
    infns = glob.glob(tlusty_dir + '/*.vis.7')
    for ii,infn in enumerate(infns):
		indata = readcol(infn)
		wav    = indata[:,0]
		data   = indata[:,1]
		cdata  = np.convolve(data,np.ones(10)/10.0,'same')
		intwav = 0.1*np.arange(90000)+3000.0
		icdata = np.interp(intwav,wav,cdata)
		n1     = infn.split('/')[-1].split('BG')[1].split('g')	
		n2     = 'g+'+str(float(n1[1].split('v')[0])/100.0)
		n1     = 'p'+str(int(n1[0])/1)
		outname = tlusty_conv_dir+'/'+n1 + ':'+n2+':m0.0:t01:z+0.00:a+0.00.TLUSTYconv.fits'
		pyfits.writeto(outname,icdata,clobber=True)
		print('convolving '+ str(ii+1) +' out of ' + str(len(infns)))

def conv_phoenix_spect(pho_dir,pho_conv_dir):    
    """
    Take all phoenix spectra from a directory, convolve to 0.1A, 
    then save to a new directory
    Currently resampling onto a wavelength grid of 0.1A also, from 
    3000 to 12000A to match AMBRE spectra
    also mostly matching filenames
    """
    infns = glob.glob(pho_dir + '/*.fits')    
    for ii,infn in enumerate(infns):
		data = pyfits.getdata(infn)
		wav  = pyfits.getdata('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
		##go from vacuum to air wavelengths
		wav = wav/(1.0+2.735182E-4+131.4182/wav**2+2.76249e8/wav**4)
		cdata  = np.convolve(data,np.ones(10)/10.0,'same')
		intwav = 0.1*np.arange(90000)+3000.0
		icdata = np.interp(intwav,wav,cdata)
		n1     = infn.split('/')[-1].split('lte')[1].split('-')
		n2     = 'g'+n1[1]
		n1     = 'p'+n1[0]
		outname = pho_conv_dir+'/'+n1 + ':'+n2+':m0.0:t01:z+0.00:a+0.00.PHOENIXconv.fits'
		pyfits.writeto(outname,icdata,clobber=True)
		print('convolving '+ str(ii+1) +' out of ' + str(len(infns)))
   
    
def make_wifes_p08_template(fn, out_dir,outfn='',rv=0.0, click=False, h_corr=None, correct_absorption=False, correct_skylines=False):
    """From a p08 file, create a template spectrum for future cross-correlation.
    The template is interpolated onto a 0.1 Angstrom grid (to match higher resolution 
    templates.
    
    Parameters
    ----------
    ddir: string
        Data directory for the p08 file
        
    fn: string
        p08 fits filename
        
    out_dir: string
        Output directory
    
    """
    flux_stamp,wave,var,sky = read_and_find_star_p08(fn, sky_rad=10,manual_click=click, return_sky=True)
    
    spectrum,sig = weighted_extract_spectrum(flux_stamp, var)
    if correct_absorption:
        abs_temps = glob.glob('/Users/rajikak/tools/absorption_spec_F.fits')
        abs_bad_ints =([0,5500],[5500,6860],[6910,7000])
        rv_abs,rv_sig_abs,temp_used = calc_rv_template(spectrum,wave,sig,abs_temps, (abs_bad_ints), heliocentric_correction=0.0)
        print("MAKING TEMPLATE AND ABSORPTION OFFSET= " + str(rv_abs) + ', '+ str(rv_sig_abs))
        if np.isnan(rv_abs) or np.isnan(rv_sig_abs):
            print("COULD NOT GET A RELIABLE ABSORPTION RV")
            rv_abs = 0.0
        if np.isnan(rv_sig_abs):
            rv_sig_abs = np.inf
    elif correct_absorption == False:
        rv_abs = 0.0

    if correct_skylines:
        sky_spectrum,sky_sig = weighted_extract_spectrum(sky, var)
        sky_intervals=([0,5500],[5700,5850],[6000,6100],[6700,7000])
        sky_temp = glob.glob('/Users/rajikak/tools/wifes_sky.fits')
        rv_sky,rv_sig_sky, temp_used = calc_rv_template(sky_spectrum,wave,sky_sig,sky_temp, (sky_intervals), heliocentric_correction=0.0)
        print("MAKING TEMPLATE AND SKY OFFSET= " + str(rv_sky) + ', '+ str(rv_sig_sky))
        if np.isnan(rv_sky) or np.isnan(rv_sig_sky) or np.abs(rv_sky)>150.:
            print("COULD NOT GET A RELIABLE SKY RV")
            rv_sky = 0.0
            rv_sig_sky = np.nan
        if np.isnan(rv_sig_sky):
            rv_sig_sky = np.inf
    elif correct_skylines == False:
        rv_sky = 0.0

    if correct_absorption and correct_skylines:
        if rv_sig_sky < rv_sig_abs:
            rv_abs = 0.0
            print("USING CORRECTION FROM SKYLINES")
        elif rv_sig_sky > rv_sig_abs:
            rv_sky = 0.0
            print("USING CORRECTION FROM ABSORPTION")

    if h_corr == None:
        heliocentric_correction = pyfits.getheader(fn)['RADVEL']
    else:
        heliocentric_correction = h_corr
    star = pyfits.getheader(fn)['OBJNAME']
    if star == 'object' or star == '':
        star = pyfits.getheader(fn)['OBJECT']
    spectrum,sig = weighted_extract_spectrum(flux_stamp, var)
    wave = wave[:-500]
    spectrum = spectrum[:-500]
    dell_template = 0.1
    wave_template=np.arange(90000)*dell_template + 3000
    spectrum_interp = np.interp(wave_template,wave*(1 - (rv - heliocentric_correction + rv_abs+rv_sky)/2.998e5),spectrum)
    if outfn == '':
        outfn = out_dir + star + '_' + fn.split('/')[-1].split('p08')[0] + 'fits'
    else:
        outfn = out_dir + outfn + '.fits'
    pyfits.writeto(outfn,spectrum_interp,clobber=True)
    print('saved template:'+outfn)
    

def rv_fit_mlnlike(shift,modft,data,errors,gaussian_offset):
    """Return minus the logarithm of the likelihood of the model fitting the data
    
    Parameters
    ----------
    shift: float
        Shift in pixels
    modft: array-like
        Real numpy Fourier transform of the model spectrum.
    data: array-like
        spectral data.
    errors: array-like
        uncertainties in spectral data
    gaussian_offset: float
        Offset to Gaussian uncertainty distribution
    """
    shifted_mod = np.fft.irfft(modft * np.exp(-2j * np.pi * np.arange(len(modft))/len(data) * shift))
    return -np.sum(np.log(np.exp(-(data - shifted_mod)**2/2.0/errors**2) + gaussian_offset))

def rv_shift_binary(shift1, shift2, alpha, modft1, modft2):
    """Shift two templates and add them, to model a binary star"""
    data_len = (len(modft1)-1)*2
    shifted_mod1 = np.fft.irfft(modft1 * np.exp(-2j * np.pi * np.arange(len(modft1))/data_len * shift1))
    shifted_mod2 = np.fft.irfft(modft2 * np.exp(-2j * np.pi * np.arange(len(modft2))/data_len * shift2))
    return (shifted_mod1 + alpha*shifted_mod2)/(1.0 + alpha)

    
def make_fake_binary(spect,wave,sig, template_fns, flux_ratio, rv0, rv1):
    """Make a fake binary in order to test todcor etc!"""
#    (wave_log, spect_int, sig_int, template_ints) =  \
#        interpolate_spectra_onto_log_grid(spect,wave,sig, template_fns)
    
    wave_templates = []
    spect_templates = []
    for template_fn in template_fns:
        dd = np.loadtxt(template_fn)
        wave_templates.append(dd[:,0])
        spect_templates.append(dd[:,1])
    wave_templates = np.array(wave_templates)
    spect_templates = np.array(spect_templates)
    
    c_light = 3e5
    fake_binary = np.interp(wave_templates[0]*(1 - rv0/c_light),wave_templates[0], spect_templates[0]) + \
                  np.interp(wave_templates[0]*(1 - rv1/c_light),wave_templates[1], spect_templates[1])*flux_ratio
    
    #fake_binary = np.interp(wave_log*(1 - rv0/c_light),wave_log, template_ints[0]) + \
    #              np.interp(wave_log*(1 - rv1/c_light),wave_log, template_ints[1])*flux_ratio
    #Un-continuum-subtract
    #binspect = fake_binary + 1
    #return binspect, wave_log, np.ones(len(binspect))*0.01
    return fake_binary, wave_templates[0], np.ones(len(wave_templates[0]))*0.01

def interpolate_spectra_onto_log_grid(spect,wave,sig, template_dir,bad_intervals=[],\
        smooth_distance=201,convolve_template=True, nwave_log=int(1e4), subtract_smoothed=True):
    """Interpolate both the target and template spectra onto a common wavelength grid"""
    
    #Create our logarithmic wavelength scale with the same min and max wavelengths as the
    #target spectrum, and nwave_log wavelengths.
    wave_log = np.min(wave)*np.exp(np.log(np.max(wave)/np.min(wave))/\
        nwave_log*np.arange(nwave_log))
    
    #Interpolate the target spectrum onto this scale
    spect_int = np.interp(wave_log,wave,spect)
    sig_int = np.interp(wave_log,wave,sig)
    
    #Normalise 
    sig_int /= np.median(spect_int)
    spect_int /= np.median(spect_int)
    
    #Remove bad intervals
    for interval in bad_intervals:
        wlo = np.where(wave_log > interval[0])[0]
        if len(wlo)==0: 
            continue
        whi = np.where(wave_log > interval[1])[0]
        if len(whi)==0:
            whi = [len(wave_log)-1]
        whi = whi[0]
        wlo = wlo[0]
        spect_int[wlo:whi] = spect_int[wlo] + np.arange(whi-wlo,dtype='float')/(whi-wlo)*(spect_int[whi] - spect_int[wlo])
        sig_int[wlo:whi]=1
    
    if subtract_smoothed:
        #Subtract smoothed spectrum
        spect_int -= spect_int[0] + np.arange(len(spect_int))/(len(spect_int)-1.0)*(spect_int[-1]-spect_int[0])
        spect_int -= np.convolve(spect_int,np.ones(smooth_distance)/smooth_distance,'same')
    
    #Now we find the interpolated template spectra, template_ints
    template_fns = template_dir
    template_ints = np.zeros( (len(template_fns),len(wave_log)) )
    for i,template_fn in enumerate(template_fns):
        try:
            #Try loading a reduced WiFeS file first... 
            if template_fn.find("p08") >=  len(template_fn) - 8:
                print('Using raw wifes p08 file')
                flux,wave_template=read_and_find_star_p08(template_fn)
                spect_template,dummy = weighted_extract_spectrum(flux)
                dell_template = np.mean(wave_template[1:]-wave_template[:-1])
            #Try loading pickled RV standards
            elif template_fn.find('pkl') >= len(template_fn)-4:
                print('Using pickled Standards')
                template_file = open(template_fn, 'r')
                wave_template, spect_template = pickle.load(template_file)
                dell_template = np.mean(wave_template[1:]-wave_template[:-1])
            #Next try a template text file (wavelength and flux in 2 columns)
            elif template_fn.find('txt') >= len(template_fn)-4:
                print('Using text file input')
                dd = np.loadtxt(template_fn)
                dell_template = np.mean(dd[1:,0]-dd[:-1,0])
                wave_template = dd[:,0]
                spect_template = dd[:,1]
            #Finally try the Ambre convolved spectral format.
            elif template_fn.find('fit') >= len(template_fn)-4:
                #print('Using ambre models (fits with fixed wavelength grid)')
                spect_template = pyfits.getdata(template_fn)
                dell_template = 0.1
                wave_template=np.arange(90000)*dell_template + 3000
            else:
                print('Invalid rv standard or model file: ' + template_fn)
                raise UserWarning
        except:
            print('Error loading model spectrum')
            raise UserWarning
            
        #Amount of subsampling in the template
        template_subsamp = int((wave[1]-wave[0])/dell_template)
        #Make sure it is an odd number to prevent shifting...
        template_subsamp = np.maximum((template_subsamp//2)*2 - 1,1)
        spect_template = np.convolve(np.convolve(spect_template,np.ones(template_subsamp)/template_subsamp,'same'),\
                                  np.ones(2*template_subsamp+1)/(2*template_subsamp+1),'same')

        #Interpolate onto the log wavelength grid.
        template_int = np.interp(wave_log,wave_template,spect_template)
        #Normalise 
        template_int /= np.median(template_int)
        #Remove bad intervals 
        for interval in bad_intervals:
            wlo = np.where(wave_log > interval[0])[0]
            if len(wlo)==0: 
                continue
            whi = np.where(wave_log > interval[1])[0]
            if len(whi)==0:
                whi = [len(wave_log)-1]
            whi = whi[0]
            wlo = wlo[0]
            template_int[wlo:whi] = template_int[wlo] + np.arange(whi-wlo, dtype='float')/(whi-wlo)*(template_int[whi] - template_int[wlo])
        if subtract_smoothed:
            #Subtract smoothed spectrum
            template_int -= template_int[0] + np.arange(len(template_int))/(len(template_int)-1.0)*(template_int[-1]-template_int[0])
            template_int -= np.convolve(template_int,np.ones(smooth_distance)/smooth_distance,'same')
        template_ints[i,:] =  template_int
        
    return wave_log, spect_int, sig_int, template_ints
    
def calc_rv_template(spect,wave,sig, template_dir,bad_intervals,smooth_distance=101, \
    gaussian_offset=1e-4,nwave_log=1e4,oversamp=1,fig_fn='',convolve_template=True,\
    starnumber=0, plotit=False, save_figures=False, save_dir='./', heliocentric_correction=0.):
    """Compute a radial velocity based on an best fitting template spectrum.
    Teff is estimated at the same time.
    
    Parameters
    ----------
    spect: array-like
        The reduced WiFeS spectrum
        
    wave: array-like
        The wavelengths corresponding to the reduced WiFeS spectrum
        
    template_conv_dir: string
        The directory containing template spectra convolved to 0.1 Angstrom resolution
        
    bad_intervals: 
        List of wavelength intervals where e.g. telluric absorption is bad.
        
    smooth_distance: float
        Distance to smooth for "continuum" correction
        
    oversamp: float
        Oversampling of the input wavelength scale. The slit is assumed 2 pixels wide.
    
    gaussian_offset: float
        Offset for the likelihood function from a Gaussian normalised to 1. 

        
    Returns
    -------
    rv: float
        Radial velocity in km/s
    rv_sig: float
        Uncertainty in radial velocity (NB assumes good model fit)
    temp: int
        Temperature of model spectrum used for cross-correlation.
    """
    if isinstance(template_dir, list):
        template_fns = template_dir
    else:
        template_fns = glob.glob(template_dir)
    
    #ADD IN HELIOCENTRIC CORRECTION SOMEWHERE:
    #Make the Heliocentric correction...
    #rv += h['RADVEL']
    
    #Interpolate the target and template spectra.
    (wave_log, spect_int, sig_int, template_ints) = interpolate_spectra_onto_log_grid(spect,wave,sig, template_fns,bad_intervals=bad_intervals, smooth_distance=smooth_distance,convolve_template=convolve_template, nwave_log=nwave_log)
        
    #Do a cross-correlation to the nearest "spectral pixel" for each template
    drv = np.log(wave_log[1]/wave_log[0])*2.998e5
    rvs = np.zeros(len(template_fns))
    peaks = np.zeros(len(template_fns))
    for i,template_fn in enumerate(template_fns):
        template_int = template_ints[i]
        if save_figures == True:
            plt.clf()
            plt.plot(wave_log, template_int, label='template')
            plt.plot(wave_log, spect_int, label='spectrum')
            plt.title('Template no.'+str(i+1))
            #plt.xlim([6800,7000])
            plt.savefig(save_dir + 'spectrum_vs_template_' + template_fns[i].split('/')[-1].split('.fits')[0] + '.eps')
            plt.clf()
        cor = np.correlate(spect_int,template_int,'same')
        ##here it's a good idea to limit where the peak Xcorrelation can be, only search for a peak within 1000 of rv=0
        ## that's and RV range of -778 to 778 for the default spacings in the code
        peaks[i] = np.max(cor[int(nwave_log/2)-100:int(nwave_log/2)+100])/np.sqrt(np.sum(np.abs(template_int)**2))
        rvs[i] = (np.argmax(cor[int(nwave_log/2)-100:int(nwave_log/2)+100])-100)*drv
        #if starnumber == 0: print('Correlating Template ' + str(i+1)+' out of ' + str(len(template_fns)))
        #if starnumber >0  : print('Correlating Template ' + str(i+1)+' out of ' + str(len(template_fns)) +' for star '+str(starnumber))
        this_rvs = drv*(np.arange(2*smooth_distance)-smooth_distance)
        correlation = cor[int(nwave_log/2)-100:int(nwave_log/2)+100]/np.sqrt(np.sum(np.abs(template_int)**2))
        best_ind = np.argmax(correlation)
        #print("best RV for template "+str(i+1)+" is "+str(this_rvs[best_ind+1] + heliocentric_correction))
        if save_figures == True:
            plt.clf()
            plt.plot(this_rvs[1:-1], correlation/np.max(correlation))
            plt.title('Correlation_with_template_'+template_fn.split('/')[-1])
            plt.show()
            plt.savefig(save_dir + 'Correlation_with_template_' + template_fn.split('/')[-1] + '.png')
            plt.clf()
    
    
    #Find the best cross-correlation.
    ix = np.argmax(peaks)
    #print("BEST TEMPLATE:"+template_fns[ix].split('/')[-1])

    #Recompute and plot the best cross-correlation
    template_int = template_ints[ix,:]
    cor = np.correlate(spect_int,template_int,'same')
    plt.clf()
    plt.plot(drv*(np.arange(2*smooth_distance)-smooth_distance), 
        cor[int(nwave_log/2)-smooth_distance:int(nwave_log/2)+smooth_distance])

    ##store the figure data for later use
    outsave = np.array([drv*(np.arange(2*smooth_distance)-smooth_distance),cor[int(nwave_log/2)-smooth_distance:int(nwave_log/2)+smooth_distance]])
    saveoutname = fig_fn.split('.png')[0] + "_figdat.pkl"
    pickle.dump(outsave,open(saveoutname,"wb"))
    
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('X Correlation')
    #plt.show()
    fn_ix = template_fns[ix].rfind('/')
    #Dodgy! Need a better way to find a name for the template.
    fn_ix_delta = template_fns[ix][fn_ix:].find(':')
    if fn_ix_delta>0:
        name = template_fns[ix][fn_ix+1:fn_ix+fn_ix_delta]
        name_string=name
        #A little messy !!!
        if name[0]=='p':
            name = name[1:]
            name_string = 'T = ' + name + ' K'
    name_string = template_fns[ix][fn_ix+1:]
    
	#pdb.set_trace()
    #Fit for a precise RV... note that minimize (rather than minimize_scalar) failed more
    #often for spectra that were not good matches.
    modft = np.fft.rfft(template_int)
    #res = op.minimize(rv_fit_mlnlike,rvs[ix]/drv,args=(modft,spect_int,sig_int,gaussian_offset))
    #x = res.x[0]
    #res = op.minimize_scalar(rv_fit_mlnlike,args=(modft,spect_int,sig_int,gaussian_offset),bounds=((rvs[ix]-1)/drv,(rvs[ix]+1)/drv))
    #x = res.x
    #fval = res.fun
    x,fval,ierr,numfunc = op.fminbound(rv_fit_mlnlike,rvs[ix]/drv-5/drv,rvs[ix]/drv+5/drv,args=(modft,spect_int,sig_int,gaussian_offset),full_output=True)
    rv = x*drv
    rv += heliocentric_correction
    ##best model 
    shifted_mod = np.fft.irfft(modft * np.exp(-2j * np.pi * np.arange(len(modft))/len(spect_int) * x))
    #pdb.set_trace()
    fplus = rv_fit_mlnlike(x+0.5,modft,spect_int,sig_int,gaussian_offset)
    fminus = rv_fit_mlnlike(x-0.5,modft,spect_int,sig_int,gaussian_offset)
    hess_inv = 0.5**2/(fplus +  fminus - 2*fval)
    if (hess_inv < 0) | (fplus < fval) | (fminus < fval):
        #If you get here, then there is a problem with the input spectrum or fitting.
        #raise UserWarning
        print("WARNING: Radial velocity fit did not work - trying again with wider range for: " + fig_fn)
        x,fval,ierr,numfunc = op.fminbound(rv_fit_mlnlike,rvs[ix]/drv-10/drv,rvs[ix]/drv+10/drv,args=(modft,spect_int,sig_int,gaussian_offset),full_output=True)
        rv = x*drv
        #print("RV ="+str(rv)+", fval ="+str(fval))
        fplus = rv_fit_mlnlike(x+0.5,modft,spect_int,sig_int,gaussian_offset)
        #print("fplus ="+str(fplus))
        fminus = rv_fit_mlnlike(x-0.5,modft,spect_int,sig_int,gaussian_offset)
        #print("fminus ="+str(fminus))
        hess_inv = 0.5**2/(fplus +  fminus - 2*fval)
        #print("hess_inv ="+str(hess_inv))
        #import pdb
        #pdb.set_trace()
        
        if (hess_inv < 0) | (fplus < fval) | (fminus < fval):
            print("WARNING: Radial velocity fit did not work, giving up with NaN uncertainty")
            print("rv =", rv)
            #rv = np.nan

    rv_sig = np.sqrt(hess_inv*nwave_log/len(spect)/oversamp)*drv
    if rv == np.nan:
        rv_sig = np.nan

    plt.title('RV, RV_sigma:' + str(rv) + ',' +str(rv_sig))
    plt.savefig(save_dir + 'Best_correlation_temp_' + template_fns[ix].split('/')[-1] + '.png')
    plt.title(name_string + ', RV = {0:4.1f}+/-{1:4.1f} km/s'.format(rv,rv_sig))
    if len(fig_fn) > 0:
        plt.savefig(fig_fn)
    plt.clf()
    plt.plot(wave_log,spect_int)
    plt.plot(wave_log,shifted_mod)
    plt.xlim([6400.0,6700.0])
    plt.title(name_string + ', RV = {0:4.1f}+/-{1:4.1f} km/s'.format(rv,rv_sig))
    if len(fig_fn) > 0:
        fig_fn_new = fig_fn.split('_xcor.png')[0] + 'fitplot.png' 
    	plt.savefig(fig_fn_new)
    #again save the figure data for use later in making nicer plots with IDL
    outsave = np.array([wave_log,spect_int,shifted_mod])
    saveoutname = fig_fn.split('_xcor.png')[0] + 'fitplot_figdat.pkl'
    pickle.dump(outsave,open(saveoutname,"wb"))
   # pdb.set_trace()
    return rv,rv_sig,template_fns[ix].split('/')[-1]
        
def calc_rv_todcor(spect,wave,sig, template_fns,bad_intervals=[],fig_fn='',\
    smooth_distance=201,convolve_template=True, alpha=0.3,\
    nwave_log=int(1e4),ncor=1000, return_fitted=False,jd=0.0,out_fn='',\
    heliocentric_correction=0, plotit=False):
    """Compute a radial velocity based on an best fitting template spectrum.
    Teff is estimated at the same time.
    
    Parameters
    ----------
    spect: array-like
        The reduced WiFeS spectrum
        
    wave: array-like
        The wavelengths corresponding to the reduced WiFeS spectrum
        
    template_fns: string
        Spectral template for star 1 and star 2 that can be read in by np.loadtxt
        
    bad_intervals: 
        List of wavelength intervals where e.g. telluric absorption is bad. For todcor,
        These can only be smoothed over.
        
    smooth_distance: float
        Distance to smooth for "continuum" correction
        
        
    Returns
    -------
    rv1: float
        Radial velocity of star 1 in km/s
    rv_sig1: float
        Uncertainty in radial velocity (NB assumes good model fit)
    rv2: float
        Radial velocity of star 1 in km/s
    rv_sig2: float
        Uncertainty in radial velocity (NB assumes good model fit)
    corpeak: float
        Correlation peak
    """  
    (wave_log, spect_int, sig_int, template_ints) =  \
        interpolate_spectra_onto_log_grid(spect,wave,sig, template_fns,\
            bad_intervals=bad_intervals, smooth_distance=smooth_distance, \
            convolve_template=convolve_template, nwave_log=nwave_log)
        
    rvs = np.zeros(len(template_fns))
    peaks = np.zeros(len(template_fns))
    drv = np.log(wave_log[1]/wave_log[0])*2.998e5
      
    #*** Next (hopefully with two templates only!) we continue and apply the TODCOR algorithm.
    
    window_width = nwave_log//20
    ramp = np.arange(1,window_width+1,dtype=float)/window_width
    window = np.ones(nwave_log)
    window[:window_width] *= ramp
    window[-window_width:] *= ramp[::-1]
    
    template_ints[0] *= window
    template_ints[1] *= window
    spect_int *= window
    
    norm1 = np.sqrt(np.sum(template_ints[0]**2))  
    norm2 = np.sqrt(np.sum(template_ints[1]**2))    
    norm_tgt = np.sqrt(np.sum(spect_int**2))
    
    #pdb.set_trace()
    c1  = np.fft.irfft(np.conj(np.fft.rfft(template_ints[0]/norm1))*np.fft.rfft(spect_int/norm_tgt))
    c1 = np.roll(c1,ncor//2)[:ncor]
    c2  = np.fft.irfft(np.conj(np.fft.rfft(template_ints[1]/norm2))*np.fft.rfft(spect_int/norm_tgt))
    c2 = np.roll(c2,ncor//2)[:ncor]
    
    #Unclear which way around this line should be. ix_c12 sign was corrected in order to 
    #give the right result with simulated data.
    c12 = np.fft.irfft(np.fft.rfft(template_ints[1]/norm2)*np.conj(np.fft.rfft(template_ints[0]/norm1)))
    c12 = np.roll(c12,ncor//2)[:ncor]
    ix = np.arange(ncor).astype(int)
    xy = np.meshgrid(ix,ix)

    #Correct the flux ratio for the RMS spectral variation. Is this needed???
    alpha_norm = alpha * norm2/norm1
    ix_c12 = np.minimum(np.maximum(xy[0]-xy[1]+ncor//2,0),ncor-1) #!!!This was the old line !!!
    #ix_c12 = np.minimum(np.maximum(xy[1]-xy[0]+ncor//2,0),ncor-1) #XXX New (temporary?) line XXX
    todcor = (c1[xy[0]] + alpha_norm*c2[xy[1]])/np.sqrt(1 + 2*alpha_norm*c12[ix_c12] + alpha_norm**2)
    
    print("Max correlation: {0:5.2f}".format(np.max(todcor)))
    #print(alpha_norm)
    #plt.plot(drv*(np.arange(nwave_log)-nwave_log//2),np.roll(c1,nwave_log//2))
    #Figure like TODCOR paper:
    #fig = plt.figure()
    #ax = fig.gca(projection='3d') 
    #ax.plot_surface(xy[0],xy[1],todcor)
    
    plt.clf()
    plt.imshow(todcor, cmap=cm.gray,interpolation='nearest',extent=[-drv*ncor/2,drv*ncor/2,-drv*ncor/2,drv*ncor/2])

    xym = np.unravel_index(np.argmax(todcor), todcor.shape)
    hw_fit = 2
    
    if (xym[0]< hw_fit) | (xym[1]< hw_fit) | (xym[0]>= ncor-hw_fit) | (xym[1]>= ncor-hw_fit):
        print("Error: TODCOR peak to close to edge!")
        raise UserWarning
    
    ix_fit = np.arange(-hw_fit, hw_fit + 1).astype(int)
    xy_fit = np.meshgrid(ix_fit,ix_fit)
    p_init = models.Gaussian2D(amplitude=np.max(todcor),x_mean=0, y_mean=0, 
        x_stddev = 50.0/drv, y_stddev = 50.0/drv)
    fit_p = fitting.LevMarLSQFitter()
    
    p = fit_p(p_init, xy_fit[0], xy_fit[1], todcor[xym[0]-hw_fit:xym[0]+hw_fit+1, 
                                                   xym[1]-hw_fit:xym[1]+hw_fit+1])

    #import pdb; pdb.set_trace()

    rv_x = drv*((p.parameters[1] + xym[1]) - ncor//2)
    rv_y = drv*((p.parameters[2] + xym[0]) - ncor//2)

    model_spect = rv_shift_binary(rv_x/drv, rv_y/drv, alpha, np.fft.rfft(template_ints[0]), np.fft.rfft(template_ints[1]))
    
    if plotit:
        (wave_log, spect_int_norm, sig_int, template_int_norm) =  \
        interpolate_spectra_onto_log_grid(spect,wave,sig, template_fns,\
            bad_intervals=bad_intervals, smooth_distance=smooth_distance, \
            convolve_template=convolve_template, nwave_log=nwave_log, \
            subtract_smoothed=False)
        model_spect_norm = rv_shift_binary(rv_x/drv, rv_y/drv, alpha, \
            np.fft.rfft(template_int_norm[0]), np.fft.rfft(template_int_norm[1]))
        model_spect_prim = rv_shift_binary(rv_x/drv, rv_y/drv, 0, \
            np.fft.rfft(template_int_norm[0]), np.fft.rfft(template_int_norm[1]))
        model_spect_sec = rv_shift_binary(rv_x/drv, rv_y/drv, 1e6, \
            np.fft.rfft(template_int_norm[0]), np.fft.rfft(template_int_norm[1]))
        ss = np.ones(5e2)/5e2
        model_ss = np.convolve(model_spect_norm, ss, mode='same')
        spect_ss = np.convolve(spect_int_norm, ss, mode='same')
        plt.clf()
        plt.plot(wave_log, model_spect_norm/model_ss, label='Joint Model')
        plt.plot(wave_log, model_spect_prim/model_ss/(1+alpha), label='Primary')
        plt.plot(wave_log, model_spect_sec/model_ss*alpha/(1+alpha), label='Secondary')
        plt.plot(wave_log, spect_int_norm/spect_ss, label='Data')
        plt.legend()
        plt.axis([3810, 5610, 0, 1.45])
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel('Flux (normalised)')
        plt.draw()
        
        #pdb.set_trace() #XXX
    
    #Compute theoretical RV uncertainties from the "Q" factors...
    errors = []
    for i,template_int in enumerate(template_ints):
        if (i==0):
            ti = template_int/(1 + alpha)
        else:
            ti = template_int*alpha/(1 + alpha)
        model_spect_deriv = (ti[1:]-ti[:-1])/(wave_log[1:]-wave_log[:-1])
        wave2_on_s = (0.5*(wave_log[1:]+wave_log[:-1]))**2/(0.5*(ti[1:]+ti[:-1]+2))
        q_factor = np.sqrt(np.mean(wave2_on_s*model_spect_deriv**2))
        photon_rv_error = 3e5/q_factor*np.median(sig_int)/np.sqrt(len(spect))
        errors.append(photon_rv_error)

    #ISSUES: 
    #1) Error (below) not computed.
    #errors = np.sqrt(np.diag(fit_p.fit_info['cov_x']))

    if len(out_fn)>0:
        outfile = open(out_fn, 'a')
        outfile.write('{0:12.4f}, {1:8.2f}, {2:8.2f}, {3:8.2f}, {4:8.2f}, {5:8.3f}\n'.
            format(jd, rv_x + heliocentric_correction, errors[0], rv_y + heliocentric_correction, errors[1], np.max(todcor)))
        outfile.close()

    if return_fitted:
        return wave_log, spect_int, model_spect
    else:
        return rv_x, errors[0], rv_y, errors[1], np.max(todcor)
    
def rv_process_dir(ddir,template_conv_dir='./ambre_conv/',standards_dir='',outfn='rvs.txt',texfn='rvs.tex',outdir='',mask_ha_emission=False):
    """Process all files in a directory for radial velocities.
    
    Parameters
    ----------
    dir: string
        Directory in which to process the WiFeS reduced spectra
    template_conf_dir: string
        Directory containing template spectra convolved to WiFeS resolution
    outfn: string
        Output filename"""
        
    if len(standards_dir)>0:
        print("WARNING: Feature not implemented yet")
        raise UserWarning
    fns = glob.glob(ddir + '/*p08.fits'  )
  	#Uncomment to test individual stars in a data-set
    #pdb.set_trace()
    #fns = fns[32:33]
    # If an out directory isn't given, use the data directory.
    if len(outdir)==0:
        outdir=ddir
    outfile = open(outdir + '/' + outfn,'w')
    outfile.write('#name,filename,ra,dec,bmsplt,mjd,rv,sig_rv,teff \n')
    texfile = open(outdir + '/' + texfn,'w')
    for iii,fn in enumerate(fns):
    	##pdb.set_trace()
        h = pyfits.getheader(fn)
        flux,wave = read_and_find_star_p08(fn,fig_fn=outdir + '/'+ h['OBJNAME'] + '.' + h['OBSID'] + '_star.png')  		
        if h['BEAMSPLT']=='RT560':
            bad_intervals = ([0,5500],[6860,7020],)
        else:
            bad_intervals = ([6862,7020],)
    ##Maybe here decide if the star is e.g. a young K/M dwarf with lots of H-alpha emission and 
        ##bad_interval out that section of spectrum this also works to remove Ae/Be emission which causes issues       
        #pdb.set_trace()
        if mask_ha_emission == True:
            simple_spec = np.sum(np.sum(flux,axis=0),axis=0)
            harng = np.where((wave > 6560.0) & (wave < 6565.0))[0] 
            crng  = np.where((wave > 6500.0) & (wave < 6520.0))[0]
            cmed  = np.median(simple_spec[crng])
            hamed = np.median(simple_spec[harng])
            scmed = np.std(simple_spec[crng])*1.253/len(crng)
            #pdb.set_trace()
            if hamed > 5.0*scmed+cmed: bad_intervals = bad_intervals+([6550,6580],)
            print('Removing H-alpha line due to emission')
    	##pdb.set_trace()
    
        spectrum,sig = weighted_extract_spectrum(flux)
        specfn = outdir + '/' + fn[fn.rfind('/')+1:] + '.spec.csv'
        specfile = open(specfn,'w')
        for i in range(len(spectrum)):
            specfile.write('{0:6.2f},{1:6.1f},{2:6.1f}\n'.format(wave[i],spectrum[i],sig[i]))
        specfile.close()
        rv,rv_sig,name = calc_rv_template(spectrum,wave,sig,template_conv_dir, bad_intervals,\
            fig_fn=outdir + '/' + h['OBJNAME'] + '.' + h['OBSID'] + '_xcor.png',starnumber=iii+1)
        outfile.write(h['OBJNAME'] + ','+fn +','+ h['RA'] + ','+ h['DEC'] + ',' + h['BEAMSPLT'] + \
         ',{0:10.3f},{1:5.1f},{2:5.1f},'.format(h['MJD-OBS'],rv,rv_sig)+name + ' \n')
        texfile.write(h['OBJNAME'] + ' & '+ h['RA'] + ' & '+ h['DEC'] + ' & ' + h['BEAMSPLT'] + \
            ' & {0:10.3f} & {1:5.1f} $\pm$ {2:5.1f} & '.format(h['MJD-OBS'],rv,rv_sig,name) + name + '\\\\ \n')
    outfile.close()
    texfile.close()
    
