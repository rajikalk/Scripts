import glob
import astropy.io.fits as pyfits
import process_stellar as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op

files = glob.glob('/Users/rajikak/Observational_Data/Reduced_Arcs/20170727/reduced_data/*p08*')

flux_0 = []
MJD = []
RV = []

#Make template of arc for each slitlet at the same pixel
print("MAKING TEMPLATES")
#read in file
temp_file = files[0]
hdu = pyfits.open(temp_file)
#Get fluxes
flux = np.array([hdu[i].data for i in range(1,13)])
#Get wavelength scale and convert to logscale
wave_0 = hdu[1].header['CRVAL1'] + np.arange(flux.shape[2])*hdu[1].header['CDELT1']
#get MJD to save
MJD_obs = hdu[0].header['MJD-OBS']
MJD.append(MJD_obs)
for slit in range(12):
    RV.append([0.0])
    flux_stamp = flux[slit]
    
    #interpolate spectrum onto log scale
    flux_0.append(flux_stamp[19])
print("MADE TEMPLATES FOR EACH SLIT")
print("NOW CALCULATING RV VARIATION FOR EACH SLIT")

for file in files[1:]:
    print("DOING FILE:", file)
    hdu = pyfits.open(file)
    flux = np.array([hdu[i].data for i in range(1,13)])
    #get wavelength scale and convert to log scale
    wave = hdu[1].header['CRVAL1'] + np.arange(flux.shape[2])*hdu[1].header['CDELT1']
    
    MJD_obs = hdu[0].header['MJD-OBS']
    MJD.append(MJD_obs)
    for slit in range(12):
        print("DOING SLIT:", slit)
        template_int = flux_0[slit]
        template_int = np.nan_to_num(template_int)
        
        flux_stamp = flux[slit][19]
        spect_int = np.nan_to_num(flux_stamp)
        import pdb
        pdb.set_trace()
        
        drv = np.log(wave_log[1]/wave_log[0])*2.998e5
        rvs = np.zeros(1)
        peaks = np.zeros(1)
        
        cor = np.correlate(spect_int,template_int,'same')
        peaks[0] = np.max(cor[int(nwave_log/2)-100:int(nwave_log/2)+100])/np.sqrt(np.sum(np.abs(template_int)**2))
        rvs[0] = (np.argmax(cor[int(nwave_log/2)-100:int(nwave_log/2)+100])-100)*drv
        this_rvs = drv*(np.arange(2*101)-101)
        correlation = cor[int(nwave_log/2)-100:int(nwave_log/2)+100]/np.sqrt(np.sum(np.abs(template_int)**2))
        
        plt.clf()
        plt.plot(drv*(np.arange(2*101)-101),
                 cor[int(nwave_log/2)-101:int(nwave_log/2)+101])
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('X Correlation')
        #plt.show()
        
        modft = np.fft.rfft(template_int)
        gaussian_offset=1e-4
        sig_int = np.ones(len(spect_int))
        x,fval,ierr,numfunc = op.fminbound(ps.rv_fit_mlnlike,rvs[0]/drv-5/drv,rvs[0]/drv+5/drv,args=(modft,spect_int,sig_int,gaussian_offset),full_output=True)
        print(x, drv, x*drv)
        rv = x*drv
        print(rv)
        shifted_mod = np.fft.irfft(modft * np.exp(-2j * np.pi * np.arange(len(modft))/len(spect_int) * x))
        
        fplus = ps.rv_fit_mlnlike(x+0.5,modft,spect_int,sig_int,gaussian_offset)
        fminus = ps.rv_fit_mlnlike(x-0.5,modft,spect_int,sig_int,gaussian_offset)
        hess_inv = 0.5**2/(fplus +  fminus - 2*fval)
        if (hess_inv < 0) | (fplus < fval) | (fminus < fval):
            #If you get here, then there is a problem with the input spectrum or fitting.
            #raise UserWarning
            #print("WARNING: Radial velocity fit did not work - trying again with wider range for: " + fig_fn)
            x,fval,ierr,numfunc = op.fminbound(ps.rv_fit_mlnlike,rvs[0]/drv-10/drv,rvs[0]/drv+10/drv,args=(modft,spect_int,sig_int,gaussian_offset),full_output=True)
            print(x, drv, x*drv)
            rv = x*drv
            print(rv)
            #print("RV ="+str(rv)+", fval ="+str(fval))
            fplus = ps.rv_fit_mlnlike(x+0.5,modft,spect_int,sig_int,gaussian_offset)
            #print("fplus ="+str(fplus))
            fminus = ps.rv_fit_mlnlike(x-0.5,modft,spect_int,sig_int,gaussian_offset)
            #print("fminus ="+str(fminus))
            hess_inv = 0.5**2/(fplus +  fminus - 2*fval)
            #print("hess_inv ="+str(hess_inv))
            #import pdb
            #pdb.set_trace()
            
            if (hess_inv < 0) | (fplus < fval) | (fminus < fval):
                print("WARNING: Radial velocity fit did not work, giving up with NaN uncertainty")
    
        rv_sig = np.sqrt(hess_inv*nwave_log/len(flux_stamp)/1)*drv
        
        RV[slit].append(rv)
        #flux_stamp = np.array([flux_stamp])
        #spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)


print("NOW THAT WE'VE GONE THROUGH ALL THE FILES, LETS PLOT THIS!")
plt.clf()
marker = ['v','^','<','>','8','s','p','*','h','D','o','x']
for slit in range(12):
    plt.scatter(MJD, RV[slit], marker=marker[slit], label=str(slit))
plt.xlabel('MJD')
plt.ylabel('RV (km/s)')
plt.ylim([-10.0, 10.0])
plt.legend(loc='best')
plt.savefig('Arc_variation.png')