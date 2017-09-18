import process_stellar as ps
import glob
import astropy.io.fits as pyfits
import csv
import numpy as np

reduced_dir = '/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/'

#Read in Standards Data
RV_standard_name = []
RV_standard_data = []

header = 0
with open('/Users/rajikak/Observational_Data/RV_standards/RV_standard.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            RV_standard_name.append(row[0])
            RV_standard_data.append(row[1:])
        else:
            header = 1

f = open('/Users/rajikak/Observational_Data/RV_standards/RV_standard_testing.csv', 'w')
f.write('Name,RV(km/s),Spectral Type,RA,DEC,file,SNR,Retrieved RV,RV_sig,Delta_RV,template\n')

images = glob.glob('/Users/rajikak/Observational_Data/RV_standards/Images/*/*')
for image in images:
    obj_name = image.split('/')[-2]
    ind = RV_standard_name.index(obj_name)
    file_dir = image.split('/')[-1].split('-')[1].split('.')[0]
    file = image.split('/')[-1].split('.png')[0]

    full_file_path = reduced_dir + file_dir + '/' + file
    
    hdu = pyfits.open(full_file_path)
    h_corr = hdu[0].header['RADVEL']

    flux,wave = ps.read_and_find_star_p08(full_file_path)
    spectrum,sig = ps.weighted_extract_spectrum(flux)
    

    template_name = '/Users/rajikak/tools/templates/' + obj_name + '_' + file.split('.p08')[0]+'.fits'

    usable_temps = glob.glob('/Users/rajikak/tools/templates/*.fits')
    if template_name in usable_temps:
        usable_temps.remove(template_name)

    rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,usable_temps, ([0,5400],[6870,6890]), save_figures=False, heliocentric_correction=h_corr)

    Delta_RV = abs(float(RV_standard_data[ind][0]) - rv)


    write_line = obj_name + ',' + ','.join(RV_standard_data[ind]) + ',' + file + ',' + str(np.mean(spectrum)) + ',' + str(rv) + ',' + str(rv_sig) + ',' + str(Delta_RV) + ',' + temp_used + '\n'
    f.write(write_line)

f.close()

