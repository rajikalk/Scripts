import process_stellar as ps
import sys
from astropy.io import fits
import os
import glob
import csv

use_self = False
click = True
reduced_dir = sys.argv[1]
fns = glob.glob(reduced_dir+'*.p08*')
f = open(reduced_dir + 'output_file.csv', 'w')
f.write('Obj_Name, RA, DEC, Obs_date, MJD, RV, RV_sigma, helio_corr, template\n')
f.close()

for fn in fns:
    f = open(reduced_dir + 'output_file.csv', 'a')
    hdu = fits.open(fn)
    Obj_name = hdu[0].header['OBJECT']
    Obs_date = hdu[0].header['DATE-OBS'].split('T')[0]
    MJD = hdu[0].header['MJD-OBS']
    h_corr = hdu[0].header['RADVEL']
    RA = hdu[0].header['RA']
    DEC = hdu[0].header['DEC']

    temp_dir = '/Users/rajikak/temp_'+Obj_name+'_'+Obs_date+'/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    #SAVE IMAGE OF STAR
    flux,wave = ps.read_and_find_star_p08(fn, fig_fn=temp_dir+Obj_name+'_'+Obs_date+'.png', manual_click=click)
    spectrum,sig = ps.weighted_extract_spectrum(flux)

    templates = glob.glob('/Users/rajikak/tools/templates/*') + glob.glob('/Users/rajikak/tools/self_templates/*')
    for temp in templates:
        if Obj_name in temp and ''.join(Obs_date.split('-')) not in temp and use_self == True:
            templates = [temp]
            print("Observed before so will compare with itself")
            break
        else:
            templates = glob.glob('/Users/rajikak/tools/templates/*')
    if len(templates) >1:
        print("Using RV standards")
    #templates = glob.glob('/Users/rajikak/tools/self_templates/*')
    #templates = glob.glob('/Users/rajikak/Observational_Data/Reduced_Wifes_Obs/20140622/*-0030*')
    '''
    for temp in templates:
        if Obj_name in temp and ''.join(Obs_date.split('-')) in temp:
            templates.remove(temp)
            print "removed template:", temp
    '''

    #PLOT TEMPLATE SPECTRA VS REAL DATA AND PLOT CROSS_CORRELATION FUNCTIONS WITH RV
    rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, ([0,5400],[6870,6890]), save_figures=True, save_dir=temp_dir, heliocentric_correction=h_corr)

    write_string = Obj_name+','+RA+','+DEC+','+str(Obs_date)+','+str(rv)+','+str(rv_sig)+','+str(h_corr)+','+str(temp_used)+'\n'
    print(write_string)
    f.write(write_string)
    f.close()
    #MAKE LATEX

f.close()