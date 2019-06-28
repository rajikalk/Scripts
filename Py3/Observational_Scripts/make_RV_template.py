import process_stellar as ps
import glob
import os
import astropy.io.fits as pyfits
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import numpy as np
import matplotlib.pyplot as plt

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-abs", "--absorption", help="do you want to correct for sky absorption?", default='False')
    parser.add_argument("-allspt", "--all_spectral_types", help="do you want to correct absorption for all spectral types", default='False')
    parser.add_argument("-sky", "--skylines", help="do you want to correct for skylines?", default='False')
    parser.add_argument("-abs_spt", "--abs_corr_spectral_types", help="Which spectral types would you like to correct the absorption line rv for? default is 'F,G'", type=str, default='F,G')
    parser.add_argument("-temp_dir", "--template_dir", help="where do you want to save the templates? if None, selects default in /tools.", default=None, type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()
if args.absorption == 'False':
    abs_cor = False
else:
    abs_cor = True
if args.skylines == 'False':
    sky_cor = False
else:
    sky_cor = True

S = Simbad()
S.add_votable_fields('rv_value')

if args.template_dir == None:
    templates_dir = '/Users/rajikak/tools/templates/'
else:
    templates_dir = args.template_dir
files = glob.glob('/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/*/*p08.fits')
image_dir = '/Users/rajikak/Observational_Data/RV_standards/Images/'
RV_standard = []
RV = []
Spec_type = []
RA = []
DEC = []
spectra_F = []
sky_spectra = []
abs_corr_spt = args.abs_corr_spectral_types.split(',')

dell_template = 0.1
wave_template=np.arange(90000)*dell_template + 3000

f = open('/Users/rajikak/Observational_Data/RV_standards/RV_standard.csv', 'w')
f.write('Name,RV(km/s),Spectral Type,RA,DEC\n')

for file in files:
    print("Doing file:", file)
    hdu = pyfits.open(file)
    Obj_name = hdu[0].header['OBJNAME']
    if Obj_name == 'object' or Obj_name == '':
        Obj_name = hdu[0].header['OBJECT']
    if os.path.exists(image_dir+Obj_name +'/'+file.split('/')[-1]+'.png'):
        if not os.path.exists(image_dir + Obj_name):
            os.makedirs(image_dir + Obj_name)

        flux,wave,var = ps.read_and_find_star_p08(file, fig_fn=image_dir+Obj_name +'/'+file.split('/')[-1]+'.png')

        #plt.clf()
        #plt.plot(wave, spectrum)
        #plt.show()
        
        if Obj_name not in RV_standard:
            RV_standard.append(Obj_name)
            
            custom_simbad = Simbad()
            custom_simbad.add_votable_fields("sptype")
            custom_simbad.add_votable_fields("rv_value")
            data = custom_simbad.query_object(Obj_name)
            
            try:
                result = Vizier.query_object(Obj_name)
                interesting_table = result['J/ApJS/141/503/table1']
                object_RV = float(interesting_table[0]['__RV_'])/1000.
            except:
                object_RV = float(data['RV_VALUE'])
            
            sptype = data['SP_TYPE'][-1]
            spectrum,sig = ps.weighted_extract_spectrum(flux,var)

            if (sptype[0] == 'F'):
                spectrum_interp = np.interp(wave_template,wave[3300:]*(1 - (0.0 - 0.0)/2.998e5),spectrum[3300:])
                spectra_F.append(spectrum_interp)
            Spec_type.append(sptype)
            RV.append(object_RV)
            ra = str(data['RA'][-1])
            RA.append(ra)
            dec = str(data['DEC'][-1])
            DEC.append(dec)
            write_string = Obj_name + ',' + str(object_RV) + ',' + sptype + ',' + ra + ',' + dec + '\n'
            f.write(write_string)

            if args.all_spectral_types != 'False':
                print("Making template with absorption correction for all spectral types")
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=abs_cor, correct_skylines=sky_cor)
            elif sptype[0] in abs_corr_spt:# or (sptype[0] == 'K'):
                print("Making template with absorption correction for selected spectral types")
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=abs_cor, correct_skylines=sky_cor)
            else:
                print("Making template with no absorption and sky=", sky_cor)
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=False, correct_skylines=sky_cor)
            
        else:
            ind = RV_standard.index(Obj_name)
            object_RV = RV[ind]
            
            if args.all_spectral_types != 'False':
                print("Making template with absorption correction for all spectral types")
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=abs_cor, correct_skylines=sky_cor)
            elif sptype[0] in abs_corr_spt:
                print("Making template with absorption correction for selected spectral types")
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=abs_cor, correct_skylines=sky_cor)
            else:
                print("Making template with no absorption and sky=", sky_cor)
                ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=False, correct_skylines=sky_cor)


        print("Object:", Obj_name, "has RV:", object_RV)

        #ps.make_wifes_p08_template(file, templates_dir, rv=object_RV, correct_absorption=True, correct_skylines=False)

f.close()

absorption_spec_F = np.median(spectra_F,axis=0)
outfn = '/Users/rajikak/tools/absorption_spec_F.fits'
pyfits.writeto(outfn,absorption_spec_F,clobber=True)

