#plot processed data
import pyfits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numpy.ma as ma
import sys
import argparse
import glob
import process_stellar as ps
import csv
import os

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="Where are the input files", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Where do you want to save output files? defaults to current dir", type=str, default='./')
    parser.add_argument("-sn", "--save_name", help="What do you want to save the output as?", type=str, default=None)
    parser.add_argument("-pltyp", "--plot_type", help="What type of plot do you want to make?", type=str, default='spec')
    parser.add_argument("-update", "--update_database", help="do you want to update the database file?", default=True)
    args = parser.parse_args()
    return args

#=============================================================================================
csv_file = '/home/rajikak/Databases/median_values.csv'
Objects = {}
args = parse_inputs()

first = True
header = True
with open(csv_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header == False:
            if len(row[0]) != 0:
                name = row[0]
            date = row[1]
            mjd = float(row[2])
            median = float(row[3])
            if (name in Objects) == False:
                Objects.update({name:np.array([[date,mjd,median]])})
            else:
                obj_data = Objects[name]
                dict_date = np.append(obj_data.T[0],date)
                dict_mjd = np.append(obj_data.T[1],mjd)
                dict_medians = np.append(obj_data.T[2],median)
                Objects[name] = np.array([dict_date, dict_mjd,dict_medians]).T
        else:
            header = False

if args.update_database:
    os.remove(csv_file)

files = sorted(glob.glob(args.input_dir + '*p08*'))

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

for fn in files:
    a = pyfits.open(fn)
    objname = a[1].header['OBJECT']
    if objname[:2] == 'UU':
        a.close()
        a = pyfits.open(fn,mode='update')
        a[1].header['OBJECT'] = objname[1:]
        a.close()
        a = pyfits.open(fn)
        objname = a[1].header['OBJECT']
    if objname.isdigit():
        a.close()
        a = pyfits.open(fn,mode='update')
        a[1].header['OBJECT'] = 'UCAC4-'+objname
        a.close()
        a = pyfits.open(fn)
        objname = a[1].header['OBJECT']
    if len(objname) == 0:
        objname = 'object'
    obsmjd = a[1].header['MJD-OBS']
    obsdate = a[1].header['DATE-OBS'].split('T')[0]
    if args.save_name == None:
        save_name = objname + '_spectra.eps'
    else:
        save_name = args.save_name
    flux,wave = ps.read_and_find_star_p08(fn)
    spectrum,sig = ps.weighted_extract_spectrum(flux)

    median_value = np.median(spectrum)
    if objname in Objects:
        obj_data = Objects[objname]
        if str(obsmjd) not in obj_data[:,1]:
            dict_date = np.append(obj_data.T[0],obsdate)
            dict_mjd = np.append(obj_data.T[1],obsmjd)
            dict_medians = np.append(obj_data.T[2],median_value)
            Objects[objname] = np.array([dict_date,dict_mjd,dict_medians]).T
    else:
        Objects.update({objname:np.array([[obsdate,obsmjd,median_value]])})

    #plot figure
    plt.clf()
    plt.plot(wave, spectrum)
    plt.xlim([np.min(wave), np.max(wave)])
    plt.xlabel('wavelength')
    plt.ylabel('flux')
    plt.savefig(save_name, bbox_inches='tight')
    print("created figure:", save_name)
    '''
    flux = np.array([a[i].data for i in range(1,13)])
    wave = a[1].header['CRVAL1'] + np.arange(flux.shape[2])*a[1].header['CDELT1']
    flux = flux.reshape(456, 6991)
    fig = plt.figure(figsize=(12, 12))
    flux_1d = flux.reshape(3187896)
    mask = ma.masked_invalid(flux_1d)
    vmax_val = np.max(flux_1d*mask)
    im = plt.imshow(flux, interpolation='nearest', vmin=0, vmax=4000, extent=[wave[0],wave[-1],456,0])
    plt.xlabel('Wavelength($\AA$)')
    plt.ylabel('row')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()
    plt.savefig(sn, bbox_inches='tight')
    '''
if args.update_database:
    f = open(csv_file, 'w')
    f.write('#Object, Date, MJD, Median \n')
    for object in range(len(list(Objects.keys()))):
        obj_name = list(Objects.keys())[object]
        obj_data = Objects[obj_name]
        for data_point in range(len(obj_data)):
            obj_date = obj_data[data_point][0]
            obj_mjd = obj_data[data_point][1]
            obj_median = obj_data[data_point][2]
            if data_point == 0:
                f.write(obj_name + ',' + str(obj_date) + ',' + str(obj_mjd) + ',' + str(obj_median) + '\n')
            else:
                f.write(',' + str(obj_date) + ',' + str(obj_mjd) + ',' + str(obj_median) + '\n')
    f.close()