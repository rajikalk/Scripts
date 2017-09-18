from astropy.io import fits
import glob
import csv

files = glob.glob('*.fits')
f = open('file_coordinates.csv', 'w')
f.write('Obj_Name, Filename, RA, DEC\n')

for file in files:
    hdu = fits.open(file)
    if hdu[0].header['imagetyp'] == 'object':
        object_name = hdu[0].header['object']
        RA = hdu[0].header['RA']
        DEC = hdu[0].header['DEC']
        f.write(object_name + ',' + file + ',' + RA + ',' + DEC + '\n')
    hdu.close()

f.close()