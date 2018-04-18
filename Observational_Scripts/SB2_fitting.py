import csv
import glob
import numpy as np
import os
from astropy.io import fits
import process_stellar as ps
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
from fpdf import FPDF
import pylab as pl
from scipy import optimize
import pickle
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import sys

def function(x, a, b, c):
    return np.exp(-(np.log(x) - b)**2 / (2 * c**2))*( a / (x * c * np.sqrt(2 * np.pi)))
            

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", "--input_file", help="what file do you want to read in?", type=str)
    parser.add_argument("-plot_curves", "--plot_curves", help="Dow you want to plot the RV curves?", type=str, default='True')
    parser.add_argument("-plot_suffix", "--plot_suffix", help="Do you want to give your plots suffixes to help distinguish between them?", type=str, default='')
    parser.add_argument("-abs_spt", "--abs_corr_spectral_types", help="Which spectral types would you like to correct the absorption line rv for? default is 'F,G,K,M'", type=str, default='G')
    parser.add_argument("-temp_dir", "--template_dir", help="what template directory do you want to use? If none, it uses the default in the /tools directory", type=str, default=None)
    args = parser.parse_args()
    return args

def main():
    
    rank = CW.Get_rank()
    size = CW.Get_size()
    
    args = parse_inputs()
    if args.template_dir == None:
        templates_dir = '/Users/rajikak/tools/templates/'
    else:
        templates_dir = args.template_dir

    #set up lists and arrays
    Object = []
    Files = []
    Other_template = []
    Other_template_SpT = []
    Alpha = []
    Pref_template = []
    Temp_sptype = []
    abs_corr_spt = args.abs_corr_spectral_types.split(',')
    Errors = []
    RVs = []
    H_corr = []

    all_templates = glob.glob(templates_dir + '/*.fits')

    if args.template_dir == None:
        templates_dir = '/Users/rajikak/tools/templates/'
    else:
        templates_dir = args.template_dir

    abs_temp = glob.glob('/Users/rajikak/tools/absorption_spec_F.fits')
    sky_temp = glob.glob('/Users/rajikak/tools/wifes_sky.fits')

    sky_intervals = ([0,5500],[5700,5850],[6000,6100],[6700,6800])
    abs_bad_ints = ([0,5500],[5500,6860],[6910,7000])

    dell_template = 0.1
    wave_template=np.arange(90000)*dell_template + 3000

    #read in current data
    print "Reading in current spreadsheet"
    header = 0
    with open('/Users/rajikak/SB2_files.csv', 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] not in Object:
                Object.append(row[0])
                append_it = len(Object)-1
                Files.append([])
                Other_template.append([])
                Other_template_SpT.append([])
                Alpha.append([])
                Errors.append([])
                RVs.append([])
                H_corr.append([])
                Pref_template.append(row[2])
                Temp_sptype.append(row[3])
            else:
                append_it = Object.index(row[0])
            Files[append_it].append(row[1])

    rit = 0
    for obj in range(len(Object)):
        for file in range(len(Files[obj])):
            if rank == rit:
                hdu = fits.open(Files[obj][file])
                Obj_name = hdu[0].header['OBJECT']
                if 'U4' in Obj_name:
                    Obj_name = 'UCAC4' + Obj_name.split('U4')[-1]
                Obs_date = hdu[0].header['DATE-OBS'].split('T')[0]
                MJD = hdu[0].header['MJD-OBS']
                h_corr = hdu[0].header['RADVEL']
                flux, wave, var, sky = ps.read_and_find_star_p08(Files[obj][file], fig_fn='/Users/rajikak/Observational_Data/Todcor/'+Obj_name+'/'+Obs_date+'_'+str(MJD)+'.png', manual_click=False, return_sky=True, sky_rad=5)
                
                pref_templates = glob.glob(templates_dir+Pref_template[obj]+'*')
                pref_temp_obj = pref_templates[0].split('_')[0].split('/')[-1]
                
                spectrum,sig = ps.weighted_extract_spectrum(flux, var)
                rv_abs, rv_abs_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,abs_temp, (abs_bad_ints), heliocentric_correction=0.0)
                temp_hdu = fits.open('/Users/rajikak/tools/'+temp_used)
                spectrum_interp = np.interp(wave_template,wave*(1 - (rv_abs)/2.998e5),spectrum)
                    
                print "Absorption offset=   " + str(rv_abs) + ',    '+ str(rv_abs_sig)
                
                if np.isnan(rv_abs) or np.isnan(rv_abs_sig) or (Temp_sptype[obj][0] not in abs_corr_spt):
                    print "Could not get reliable absorption RV or wrong spectral type"
                    rv_abs = 0.0
                    rv_abs_sig = np.nan
                
                if np.isnan(rv_abs_sig):
                    rv_abs_sig = np.inf
                
                rv_sky = 0.0
                rv_err = rv_abs_sig
                print("USING CORRECTION FROM ABSORPTION:    " + str(rv_abs) + ",    " + str(rv_abs_sig))
                
                rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,pref_templates, ([0,5400],[6550, 6600],[6870,6890]), save_figures=False, save_dir='/Users/rajikak/Observational_Data/Todcor/'+Obj_name+'/', heliocentric_correction=(h_corr-rv_abs-rv_sky))
                
                if np.isinf(rv_err) == False:
                    rv_sig = np.sqrt(np.square(rv_sig) + np.square(rv_err))
                print "***************************************************"
                print "*"
                print "*"
                print "RV for object:", Obj_name, "is", rv, rv_sig
                print "*"
                print "*"
                print "***************************************************"
                
                errors = [0,0]
                rvs = [0,0]
                save_h_corr = h_corr-rv_abs-rv_sky
                spectrum,sig = ps.weighted_extract_spectrum(flux, var)
                best_fit = 0
                best_alpha = 0
                other_temp = ''
                other_temp_sptype = ''
                for template in all_templates:
                    temp_name = template.split('_')[0].split('/')[-1]
                    custom_simbad = Simbad()
                    custom_simbad.add_votable_fields("sptype")
                    data = custom_simbad.query_object(temp_name)
                    sptype = data['SP_TYPE'][-1]
                    #if pref_temp_obj != temp_name and Temp_sptype[obj] != sptype:
                    for alpha in np.linspace(0,1,101)[1:]:
                        try:
                            rv_x, errors_x, rv_y, errors_y, peak = ps.calc_rv_todcor(spectrum, wave,sig, [pref_templates[0], template], bad_intervals=[[0,5400],[6550, 6600],[6870,6890]], alpha=alpha, plotit=False)
                        except:
                            rv_x = 0
                            errors_x = 0
                            rv_y = 0
                            errors_y = 0
                            peak = 0
                        if peak > best_fit:
                            print "Current best fit:", best_fit, "with templates:", pref_templates[0], other_temp
                            best_fit = peak
                            best_alpha = alpha
                            save_h_corr = save_h_corr
                            errors = [errors_x, errors_y]
                            rvs = [rv_x, rv_y]
                            other_temp = template
                            other_temp_sptype = sptype

                send_data = [obj,other_temp.split('/')[-1],other_temp_sptype,best_alpha,errors,rvs,save_h_corr]
                if rank == 0:
                    fitting_update = send_data
                    Other_template[fitting_update[0]].append(fitting_update[1])
                    Other_template_SpT[fitting_update[0]].append(fitting_update[2])
                    Alpha[fitting_update[0]].append(fitting_update[3])
                    Errors[fitting_update[0]].append(fitting_update[4])
                    RVs[fitting_update[0]].append(fitting_update[5])
                    H_corr[fitting_update[0]].append(fitting_update[6])
                    f = open('best_fitting_file.csv', 'a')
                    write_string = Object[fitting_update[0]] + ',' + Pref_template[obj] + ',' + fitting_update[1] + ',' + Temp_sptype[obj] + ',' + fitting_update[2] + ',' + str(fitting_update[3]) + ',' + str(fitting_update[4]) + ',' + str(fitting_update[5]) + ',' + str(fitting_update[6]) '\n'
                    f.write(write_string)
                    f.close()
                    del fitting_update
                    del write_string
                else:
                    CW.send(send_data, dest=0, tag=rank)
                print "***************************************************"
                print "For object:", Object[obj]
                print "Best Alpha:", best_alpha
                print "RVs:", rvs
                print "RV errors:", errors
                print "Heliocentric Correction:", h_corr-rv_abs-rv_sky
                print "Preferred Template:", pref_templates[0], "of SpT:", Temp_sptype[obj]
                print "Other Template:", other_temp, "of SpT:", other_temp_sptype
                print "***************************************************"

                rv_x, errors_x, rv_y, errors_y, peak = ps.calc_rv_todcor(spectrum, wave,sig, [pref_templates[0], other_temp], bad_intervals=[[0,5400],[6550, 6600],[6870,6890]], alpha=best_alpha, plotit=True, fig_fn='/Users/rajikak/Observational_Data/Todcor/'+Obj_name+'/'+str(MJD)+'.png')

            rit = rit + 1
            if rit == size:
                sys.stdout.flush()
                CW.Barrier()
                rit = 0
                if rank == 0:
                    print "UPDATING CALCULATED BAYES VALUES"
                    for orit in range(1,size):
                        fitting_update = CW.recv(source=orit, tag=orit)
                        Other_template[fitting_update[0]].append(fitting_update[1])
                        Other_template_SpT[fitting_update[0]].append(fitting_update[2])
                        Alpha[fitting_update[0]].append(fitting_update[3])
                        Errors[fitting_update[0]].append(fitting_update[4])
                        RVs[fitting_update[0]].append(fitting_update[5])
                        H_corr[fitting_update[0]].append(fitting_update[6])
                        print "Updated Bayes factors retrieved from rank", orit, "for object", Object[fitting_update[0]]
                        f = open('best_fitting_file.csv', 'a')
                        write_string = Object[fitting_update[0]] + ',' + Pref_template[obj] + ',' + fitting_update[1] + ',' + Temp_sptype[obj] + ',' + fitting_update[2] + ',' + str(fitting_update[3]) + ',' + str(fitting_update[4]) + ',' + str(fitting_update[5]) + ',' + str(fitting_update[6]) '\n'
                        f.write(write_string)
                        f.close()
                        del fitting_update
                        del write_string
                sys.stdout.flush()
                CW.Barrier()

if __name__ == '__main__': main()

