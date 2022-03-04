import process_stellar as ps
import glob
import astropy.io.fits as pyfits
import csv
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

reduced_dir = '/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/'
files = glob.glob('/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/*/*p08*')
templates = glob.glob('/Users/rajikak/tools/templates/*.fits')

RV_standard = []
Actual_RV = []
Temp_to_use = []
Time_line = []
RV_variation = []

for file in files:
    print('DOING FILE:', file)
    hdu = pyfits.open(file)
    Obj_name = hdu[0].header['OBJNAME']
    if Obj_name == 'object' or Obj_name == '':
        Obj_name = hdu[0].header['OBJECT']
    MJD = hdu[0].header['MJD-OBS']

    if Obj_name not in RV_standard:
        RV_standard.append(Obj_name)
        Time_line.append([MJD])
        RV_variation.append([0.0])
        Temp_to_use.append([item for item in templates if Obj_name in item][0])
        
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
        Actual_RV.append(object_RV)
    else:
        h_corr = hdu[0].header['RADVEL']
        ind = RV_standard.index(Obj_name)
        pref_temp = Temp_to_use[ind]

        flux,wave = ps.read_and_find_star_p08(file)
        spectrum,sig = ps.weighted_extract_spectrum(flux)

        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,[pref_temp], ([0,5400],[6870,6890]), save_figures=False, heliocentric_correction=h_corr)

        Delta_RV = Actual_RV[ind] - rv

        Time_line[ind].append(MJD)
        RV_variation[ind].append(Delta_RV)

for sit in range(len(RV_standard)):
    plt.clf()
    plt.plot(Time_line[sit], RV_variation[sit], 'o')
    plt.axhline(y=0.0)
    plt.xlabel('MJD')
    plt.ylabel('RV Variation')
    plt.title(RV_standard[sit] + '_RV_' + str(Actual_RV[sit]))
    plt.savefig('/Users/rajikak/Observational_Data/RV_standards/Images/'+RV_standard[sit]+'/RV_variation.png')
