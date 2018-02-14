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

def function(x, a, b, c):
    return np.exp(-(np.log(x) - b)**2 / (2 * c**2))*( a / (x * c * np.sqrt(2 * np.pi)))
            

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", "--input_file", help="what file do you want to read in?", type=str)
    parser.add_argument("-m", "--mode", help="do you want to update the spreadsheet, or just read it in", default="update")
    parser.add_argument("-click", "--click", help="do you want to manually click?", default=False)
    parser.add_argument("-obj", "--obj_name", help="If you only want to do one object, supply name", type=str, default="")
    parser.add_argument("-start_date", "--start_d", help="what date do you want to start updating from", type=str, default='')
    parser.add_argument("-end_date", "--end_d", help="What date do you want to finish updating at?", type=str, default='')
    parser.add_argument("-date", "--date", help="What specific date do you want to do? format: YYYY-MM-DD", type=str, default='')
    parser.add_argument("-cor_sky", "--correct_sky", help="Do you want to correct for the sky?", type=str, default='True')
    parser.add_argument("-abs_only", "--absorption", help="Correct for absorption only?", type=str, default='False')
    parser.add_argument("-sky_only", "--skylines", help="Correct for sky emission only?", type=str, default='False')
    parser.add_argument("-plot_curves", "--plot_curves", help="Dow you want to plot the RV curves?", type=str, default='True')
    parser.add_argument("-plot_suffix", "--plot_suffix", help="Do you want to give your plots suffixes to help distinguish between them?", type=str, default='')
    parser.add_argument("-abs_spt", "--abs_corr_spectral_types", help="Which spectral types would you like to correct the absorption line rv for? default is 'F,G,K,M'", type=str, default='F,G')
    parser.add_argument("-temp_dir", "--template_dir", help="what template directory do you want to use? If none, it uses the default in the /tools directory", type=str, default=None)
    parser.add_argument("-hist", "--hist_pickle", help="Do you have the RV distribution already pickled?", type=str, default=None)
    args = parser.parse_args()
    return args

args = parse_inputs()

#set up lists and arrays
Object = []
Region = []
Coords = []
Membership_probability = []
IR_excess = []
Disk_visual = []
Magnitudes = []
No_obs = []
SB1_flag = []
SB2_flag = []
Mean_RV = []
RV_variation = []
Temp_sptype = []
Pref_template = []
Obs_info = []
SNR = []
abs_corr_spt = args.abs_corr_spectral_types.split(',')

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
with open(args.input_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if 'U4' in row[0]:
                row[0] = 'UCAC4' + row[0].split('U4')[-1]
            Object.append(row[0])
            if not os.path.exists('/Users/rajikak/Observational_Data/PDF_dirs/' + row[0]):
                os.makedirs('/Users/rajikak/Observational_Data/PDF_dirs/' + row[0])
            Region.append(row[1])
            Coords.append((row[2], row[3]))
            Membership_probability.append(int(row[4]))
            IR_excess.append(row[5])
            Disk_visual.append(row[6])
            Magnitudes.append([float(row[7]), float(row[8]), float(row[9])])
            No_obs.append(int(row[10]))
            if row[11] == 'FALSE' or row[11] == 'False':
                SB1_flag.append(False)
            else:
                SB1_flag.append(True)
            if row[12] == 'FALSE' or row[12] == 'False':
                SB2_flag.append(False)
            else:
                SB2_flag.append(True)
            Mean_RV.append(row[13])
            RV_variation.append(row[14])
            Pref_template.append(row[15])
            Temp_sptype.append(row[16])
            if len(row) > 17:
                Obs = np.array(row[17:])
                Obs = np.delete(Obs, np.where(Obs==''))
                #if len(Obs) > 5:
                #    Obs = np.reshape(Obs, (len(Obs)/5, 5))
                Obs = np.reshape(Obs, (len(Obs)/6, 6))
                for ind_obs in Obs:
                    if '/' in ind_obs[0]:
                        new_format = '20' + ind_obs[0].split('/')[-1] + '-' + ind_obs[0].split('/')[-2] + '-' + ("%02d" % int(ind_obs[0].split('/')[-3]))
                        ind_obs[0] = new_format
            else:
                Obs = np.array([])
            Obs_info.append(Obs)
        if header == 0:
            header = 1

Manual_click_flag = []
with open('/Users/rajikak/Observational_Data/manual_click_details.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        obj = row[0]
        date = row[1]
        selection_note = row[2]
        if row[3] != '':
            maxpx = [int(row[3]),int(row[4])]
        else:
            maxpx=[]
        Manual_click_flag.append([obj, date, selection_note, maxpx])

maxpx_list = [[],[]]
with open('/Users/rajikak/Observational_Data/manual_click_details_new_format.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        fn = row[0]
        maxpx_list[0].append(fn)
        coords = [int(row[-2]),int(row[-1])]
        maxpx_list[1].append(coords)

ignore_obs_list = []
with open('/Users/rajikak/Observational_Data/ignore_obs.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        ignore_obs_list.append([row[0], row[1]])

#Manual_click_flag = [('RIK-28', '2015-04-05', 'The Obvious Star'), ('RIK-65', '2015-05-02', 'The star closest to the center'), ('RIK-65', '2016-04-20', 'The star on the right'), ('RIK-82', '2016-03-25', 'The vague star in the middle'), ('RIK-96', '2015-05-03', 'The vague star in the middle'), ('RIK-96', '2017-07-08', 'The star on the edge'), ('RIK-146', '2017-07-08', 'The star on the edge'), ('UCAC4-4552391', '2016-03-26', 'The vague star in the middle'), ('UCAC4-39740614', '2017-07-08', 'The star on the edge'), ('UCAC4-60698038', '2017-07-06', 'The star in the middle'), ('UCAC4-60698038', '2017-07-27', 'The star in the middle'), ('UCAC4-371337294', '2016-04-18', 'The fainter smear'), ('UCAC4-379570976', '2017-07-11', 'The fainter star'), ('UCAC4-379570976', '2017-07-26', 'The fainter star'), ('UCAC4-426647186', '2017-06-02', 'The fainter of the two middle stars'), ('UCAC4-426647186', '2017-07-16', 'The fainter of the two middle stars'), ('UCAC4-44527823', '2016-04-18', 'The  star on the left of the middle two stars'), ('UCAC4-456216791', '2016-06-24', 'The faint star in the middle'), ('UCAC4-456216791', '2016-08-01', 'The faint star on the right'), ('UCAC4-469109600', '2016-03-23', 'The vague star in the middle'), ('UCAC4-817319395', '2016-04-18', 'The bright star on the top right'), ('UCAC4-827387234', '2017-06-02', 'The fainter star on the bottom left'), ('UCAC4-827387234', '2017-06-04', 'The fainter star on the bottom left'), ('UCAC4-445278523', '2016-04-18', 'Left one of the middle stars')]

#save_plots = ['RIK-96', 'UCAC4-161328427']
save_plots = []

#Reduced dirs
if args.mode == 'update' or (args.mode != 'update' and len(save_plots) > 0):
    p08_dirs = glob.glob('/Users/rajikak/Observational_Data/Reduced_Wifes_Obs/2*/')
    if args.start_d != '':
        start_ind = p08_dirs.index('/Users/rajikak/Observational_Data/Reduced_Wifes_Obs/'+args.start_d+'/')
    else:
        start_ind = 0
    if args.end_d != '':
        end_ind = p08_dirs.index('/Users/rajikak/Observational_Data/Reduced_Wifes_Obs/'+args.end_d+'/')
    else:
        end_ind = len(p08_dirs)
    p08_dirs = p08_dirs[start_ind:end_ind]
    #getting entire baseline
    start_hdu = fits.open(glob.glob(p08_dirs[0]+'*.p08*')[0])
    start_time = start_hdu[0].header['MJD-OBS']
    end_hdu = fits.open(glob.glob(p08_dirs[-1]+'*.p08*')[-1])
    end_time = end_hdu[0].header['MJD-OBS']
    for dir in p08_dirs:
        fns = glob.glob(dir+'*.p08*')
        for fn in fns:
            hdu = fits.open(fn)
            Obj_name = hdu[0].header['OBJECT']
            if 'U4' in Obj_name:
                Obj_name = 'UCAC4' + Obj_name.split('U4')[-1]
            if Obj_name == 'RIK-96' or Obj_name == 'UCAC4-161328427' or Obj_name == 'UCAC4-1253626396' or Obj_name == 'UCAC4-447414452' or Obj_name == 'UCAC4-450968247':
                ind = Object.index(Obj_name)
                f = open('SB2_files.csv', 'a')
                f.write(Obj_name + ',' + fn + ',' + Pref_template[ind] + ',' + Temp_sptype[ind] + '\n')
                f.close()
            Obs_date = hdu[0].header['DATE-OBS'].split('T')[0]
            MJD = hdu[0].header['MJD-OBS']
            #alt_date = Obs_date.split('-')[-1]+'/'+Obs_date.split('-')[-2]+'/'+Obs_date.split('-')[-3][2:]
            
            if [Obj_name, Obs_date] in ignore_obs_list:
                print "SKIPPING FILE:", fn, " BECAUSE DATA IS NOT GOOD"
            elif (Obj_name in Object and args.obj_name == "") or (Obj_name == args.obj_name) or ([item for item in Manual_click_flag if item[0] == Obj_name and item[1] == '/'.join([str(int(Obs_date.split('-')[2])), Obs_date.split('-')[1], Obs_date.split('-')[0][2:]])] != []) or (Obj_name in save_plots):
                ind = Object.index(Obj_name)
                if (str(np.round(MJD, decimals=5)) not in Obs_info[ind]) or (Obs_date == args.date) or ((Obj_name, Obs_date) in Manual_click_flag) or (Obj_name in save_plots):
                    h_corr = hdu[0].header['RADVEL']
                    RA = hdu[0].header['RA']
                    DEC = hdu[0].header['DEC']
                    
                    #determining whether you need to click on the star.
                    if args.click == 'True' or [item for item in Manual_click_flag if item[0] == Obj_name and item[1] == '/'.join([str(int(Obs_date.split('-')[2])), Obs_date.split('-')[1], Obs_date.split('-')[0][2:]])] != []:
                        click = True
                        manual_click_data = [item for item in Manual_click_flag if item[0] == Obj_name and item[1] == '/'.join([str(int(Obs_date.split('-')[2])), Obs_date.split('-')[1], Obs_date.split('-')[0][2:]])]
                        if manual_click_data[-1][-1] != []:
                            click = False
                        print "SELECTION NOTE:", manual_click_data[0][-2]
                    else:
                        click = False
                        manual_click_data = [[[]]]
                    
                    #Correct wavelenght scale using skylines if possible:
                    if manual_click_data[-1][-1] != [] or fn in maxpx_list[0]:
                        if fn in maxpx_list[0]:
                            max_px = maxpx_list[1][maxpx_list[0].index(fn)]
                        else:
                            max_px = manual_click_data[-1][-1]
                        flux, wave, var, sky = ps.read_and_find_star_p08(fn, fig_fn='/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/'+Obs_date+'_'+str(MJD)+'.png', manual_click=click, return_sky=True, sky_rad=5, max_px=max_px)
                    elif click == True:
                        flux, wave, var, sky, maxpx = ps.read_and_find_star_p08(fn, fig_fn='/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/'+Obs_date+'_'+str(MJD)+'.png', manual_click=click, return_sky=True, sky_rad=5, return_maxpx=True)
                        click_ind = Manual_click_flag.index(manual_click_data[0])
                        Manual_click_flag[click_ind] = [fn, manual_click_data[0][-2], maxpx]
                    else:
                        flux, wave, var, sky = ps.read_and_find_star_p08(fn, fig_fn='/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/'+Obs_date+'_'+str(MJD)+'.png', manual_click=click, return_sky=True, sky_rad=5)
                    
                    #selecting template
                    if Pref_template[ind] != '':
                        templates = glob.glob(templates_dir+Pref_template[ind]+'*')
                        if Temp_sptype[ind] == '':
                            temp_name = Pref_template[ind].split('_')[0]
                            custom_simbad = Simbad()
                            custom_simbad.add_votable_fields("sptype")
                            data = custom_simbad.query_object(temp_name)
                            sptype = data['SP_TYPE'][-1]
                            Temp_sptype[ind] = sptype
                        else:
                            sptype = Temp_sptype[ind]
                    else:
                        templates = glob.glob(templates_dir + '/*.fits')
                        Temp_sptype[ind] = ' '
                    
                    if args.correct_sky == 'True':
                        #Calculating RV from skylines
                        if args.absorption == 'False':
                            spectrum,sig = ps.weighted_extract_spectrum(sky,var)
                            snr_val = np.max(spectrum)/np.median(spectrum)
                            SNR.append(snr_val)
                            rv_sky,rv_sky_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,sky_temp, (sky_intervals), heliocentric_correction=0.0)
                            temp_hdu = fits.open('/Users/rajikak/tools/'+temp_used)
                            spectrum_interp = np.interp(wave_template,wave*(1 - (rv_sky)/2.998e5),spectrum)
                            
                            print "Skyline offset=  " + str(rv_sky) + ',    '+ str(rv_sky_sig)
                            
                            if np.isnan(rv_sky) or np.isnan(rv_sky_sig) or np.abs(rv_sky) > 50.:
                                print "Could not get reliable sky RV"
                                rv_sky = 0.0
                                rv_sky_sig = np.inf
                            
                            if np.isnan(rv_sky_sig) or np.abs(rv_sky) > 100.0:
                                rv_sky_sig = np.inf
                        else:
                            rv_sky = 0.0
                            rv_sky_sig = 0.0
                            which_sky = 'abs'
                    
                        if args.skylines == 'False':
                            #Calculating RV from absorption lines
                            spectrum,sig = ps.weighted_extract_spectrum(flux, var)
                            rv_abs, rv_abs_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,abs_temp, (abs_bad_ints), heliocentric_correction=0.0)
                            temp_hdu = fits.open('/Users/rajikak/tools/'+temp_used)
                            spectrum_interp = np.interp(wave_template,wave*(1 - (rv_abs)/2.998e5),spectrum)
                            
                            print "Absorption offset=   " + str(rv_abs) + ',    '+ str(rv_abs_sig)
                       
                            if np.isnan(rv_abs) or np.isnan(rv_abs_sig) or (Temp_sptype[ind][0] not in abs_corr_spt):
                                print "Could not get reliable absorption RV or wrong spectral type"
                                rv_abs = 0.0
                                rv_abs_sig = np.nan
                            
                            if np.isnan(rv_abs_sig):
                                rv_abs_sig = np.inf
                        else:
                            rv_abs = 0.0
                            rv_abs_sig = 0.0
                            which_sky = 'sky'
                    else:
                        spectrum,sig = ps.weighted_extract_spectrum(flux, var)
                        rv_sky = 0.0
                        rv_abs = 0.0
                        rv_sky_sig = 0.0
                        rv_abs_sig = 0.0
                        rv_err = 0.0
                        which_sky = ''

                    #Selecting
                    if which_sky == 'sky': #rv_sky_sig < rv_abs_sig or :
                        rv_abs = 0.0
                        rv_err = rv_sky_sig
                        print("USING CORRECTION FROM SKYLINES:  " + str(rv_sky) + ",    " + str(rv_sky_sig))
                    elif which_sky == 'abs':
                        rv_sky = 0.0
                        rv_err = rv_abs_sig
                        print("USING CORRECTION FROM ABSORPTION:    " + str(rv_abs) + ",    " + str(rv_abs_sig))

                    spectrum,sig = ps.weighted_extract_spectrum(flux, var)
                    plt.clf()
                    plt.plot(wave, spectrum)
                    plt.xlabel('wavelength')
                    plt.ylabel('flux')
                    plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/Spectrum_'+str(MJD)+'.png')
                    if Obj_name in save_plots:
                        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, ([0,5400],[6550, 6600],[6870,6890]), save_figures=True, save_dir='/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/', heliocentric_correction=(h_corr-rv_abs-rv_sky))
                        corr_pickle = '/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/cross_correlation.pkl'
                        os.rename(corr_pickle, corr_pickle.split('.pkl')[0] + '_' + str(MJD) + '.pkl')
                        
                    else:
                        rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, ([0,5400],[6550, 6600],[6870,6890]), save_figures=False, save_dir='/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/', heliocentric_correction=(h_corr-rv_abs-rv_sky))
                    
                    if Pref_template[ind] == '':
                        temp_name = temp_used.split('_')[0]
                        custom_simbad = Simbad()
                        custom_simbad.add_votable_fields("sptype")
                        data = custom_simbad.query_object(temp_name)
                        sptype = data['SP_TYPE'][-1]
                        Temp_sptype[ind] = sptype
                
                    if np.isinf(rv_err) == False:
                        rv_sig = np.sqrt(np.square(rv_sig) + np.square(rv_err))
                    print "***************************************************"
                    print "*"
                    print "*"
                    print "RV for object:", Obj_name, "is", rv, rv_sig
                    print "*"
                    print "*"
                    print "***************************************************"
                    
                    if args.plot_curves == 'True':
                        corr_img = glob.glob('/Users/rajikak/Observational_Data/PDF_dirs/'+Obj_name+'/Best_correlation_temp_*')
                        os.rename(corr_img[0], (Obs_date + '_' + str(MJD)).join(corr_img[0].split('Best')))
                    
                    if (args.date != '') or (str(np.round(MJD, decimals=5)) in Obs_info[ind]):
                        update_ind = np.where(Obs_info[ind] == str(np.round(MJD, decimals=5)))[0][0]
                        Obs_info[ind][update_ind] = [Obs_date, str(np.round(MJD, decimals=5)), rv, rv_sig, rv_abs+rv_sky, temp_used]
                    else:
                        Obs_info[ind] = np.append(Obs_info[ind],[Obs_date, str(np.round(MJD, decimals=5)), rv, rv_sig, rv_abs+rv_sky, temp_used])
                        Obs_info[ind] = np.reshape(Obs_info[ind], (len(Obs_info[ind])/6, 6))
                    No_obs[ind] = len(Obs_info[ind])

            elif Obj_name not in Object:
                print "SKIPPING FILE:", fn, " WITH OBJECT", Obj_name, "CHECK IF TYPO OR OTHER OBJECT"

f = open('/Users/rajikak/Observational_Data/manual_click_details.csv', 'w')

for obj in Manual_click_flag:
    if str(obj[-1]) == '[]':
        write_line = ','.join(obj[:-1]) + ',,\n'
    else:
        write_line = ','.join(obj[:-1]) + ',' + str(obj[-1][0]) + ',' + str(obj[-1][1]) + '\n'
    f.write(write_line)

f.close()


#Make RV plots
RV_dist = [[],[]]
RV_dist_sptype = [[], [], [], []]
#RV_dist = []
if args.hist_pickle != None:
    RV_dist = pickle.load(open(args.hist_pickle, "rb" ))
if args.mode == 'update':
    for obj in range(len(Object)):
        if No_obs[obj] != 0:
            mean_rv = np.average(Obs_info[obj][:,2].astype(np.float))
        else:
            mean_rv = 0.0
        Mean_RV[obj] = mean_rv
        if No_obs[obj] > 1:
            if args.mode == "update":
                extrema_ints = [np.argmax(Obs_info[obj][:,2].astype(np.float)), np.argmin(Obs_info[obj][:,2].astype(np.float))]
                errorbar_extrema = [float(Obs_info[obj][:,3][extrema_ints[0]]), float(Obs_info[obj][:,3][extrema_ints[1]])]
                extrema_values = [float(Obs_info[obj][:,2][extrema_ints[0]]), float(Obs_info[obj][:,2][extrema_ints[1]])]
                spread_val = extrema_values[0] - extrema_values[1]
                combined_error = np.sqrt(np.square(errorbar_extrema[0]) + np.square(errorbar_extrema[1]))
                if spread_val > 5*combined_error:
                    SB1_flag[obj] = True
                else:
                    SB1_flag[obj] = False
            else:
                spread_val = RV_variation[obj]
        
            if (Region[obj]) == 'US' and ('Y' in IR_excess[obj]):
                RV_dist[0].append(spread_val)
            elif (Region[obj]) == 'UCL' and ('Y' in IR_excess[obj]):
                RV_dist[1].append(spread_val)

            print "RV DIST FOR OBJ:", Object[obj], "IS:", spread_val
            
            if Temp_sptype[obj][0] == 'F':
                RV_dist_sptype[0].append(spread_val)
            elif Temp_sptype[obj][0] == 'G':
                RV_dist_sptype[1].append(spread_val)
            elif Temp_sptype[obj][0] == 'K':
                RV_dist_sptype[2].append(spread_val)
            elif Temp_sptype[obj][0] == 'M':
                RV_dist_sptype[3].append(spread_val)
            
            t = Obs_info[obj][:,1].astype(np.float)
            data = Obs_info[obj][:,2].astype(float)
            t_more = np.linspace(np.min(t), np.max(t), 10000)
            
            if args.plot_curves == 'True':
                plt.clf()
                err = Obs_info[obj][:,3].astype(float)
                plt.errorbar(t, data, yerr=err, fmt='o')
                '''
                guess_mean = (np.max(data)+np.min(data))/2.
                guess_std = (np.max(data)-np.min(data))/2.
                guess_phase = 0
                guess_freq = 4.*(np.max(t) - np.min(t))
                
                if len(data) > 3.:
                    data_first_guess = guess_std*np.sin(((2.*np.pi*(t_more-t_more[0]))/guess_freq)+guess_phase) + guess_mean
                    plt.plot(t_more, data_first_guess, label='first guess')
                    
                    optimize_func = lambda x: x[0]*np.sin(((2.*np.pi*t)/x[1])+x[2]) + x[3] - data
                    
                    est_std, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_std,guess_freq, guess_phase, guess_mean])[0]
                    data_fit = est_std*np.sin((2.*np.pi*est_freq*t_more)+est_phase) + est_mean
                    plt.plot(t_more, data_fit, label='after fitting')
                '''
                plt.title(Object[obj] + '_' + str(Coords[obj]))
                plt.legend(loc='best')
                plt.xlabel('MJD')
                plt.ylabel('Radial Velocity (km/s)')
                plt.xlim([start_time, end_time])
                #plt.ylim([-20.0, 20.0])
                plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/'+Object[obj]+'/RV_vs_MJD_'+Pref_template[obj]+'.png')
                print "CREATED RV CURVE FOR", Object[obj]

            RV_variation[obj] = spread_val
        else:
            RV_variation[obj] = 0.0

        #Create compiled pdf
        '''
        images = glob.glob('/Users/rajikak/Observational_Data/PDF_dirs/'+Object[obj]+'/*')
        pdf = FPDF()
        pdf.add_page()
        
        x,y,w,h,row,column = 10,10,100,68,0,0
        trim = 0.2
        for image in images:
            if '.pdf' not in image:
                if 'RV_vs_MJD' in image:
                    x = 0
                    w = w*2
                    h = h*2
                if (y+h) > 300:
                    pdf.add_page()
                    x,y,row,column = 10,10,0,0
                if x == 10 and image != images[-1]:
                    pdf.image(image, x, (y+(trim*h)), w, (h-(trim*2.*h)))
                else:
                    pdf.image(image, x, y, w, h)
                column = column + 1
                x = x + w
                if column == 2:
                    column = 0
                    x = 10
                    y = y + h

        pdf.output('/Users/rajikak/Observational_Data/PDF_dirs/'+Object[obj]+'/'+Object[obj]+'_all.pdf', 'F')
        pdf.close()
        '''

#write out updated data
print "REWRITING SPREADSHEET"
f = open(args.input_file, 'w')
f.write('Object,Region,RA,DEC,Pmem,Disk,Visual inspection,K_mag,V_mag,B_mag,No. Obs,SB1_flag,SB2 flag,Mean_RV,RV_variation,Pref_temp,Temp_SpT,Observations,MJD,RV,RV_err,RV_sky,Template\n')

for obj in range(len(Object)):
    write_string = Object[obj] + ',' + Region[obj] + ',' + Coords[obj][0] + ',' + Coords[obj][1] + ',' + str(Membership_probability[obj]) + ',' + IR_excess[obj] + ',' + Disk_visual[obj] + ',' + str(Magnitudes[obj][0]) + ',' + str(Magnitudes[obj][1]) + ',' + str(Magnitudes[obj][2]) + ',' + str(No_obs[obj]) + ',' + str(SB1_flag[obj]) + ','+ str(SB2_flag[obj]) + ',' + str(Mean_RV[obj]) + ',' + str(RV_variation[obj]) + ',' + Pref_template[obj] + ',' + Temp_sptype[obj] + ',' + ','.join(np.ravel(Obs_info[obj])) + '\n'
    f.write(write_string)

f.close()

if RV_dist[0] == []:
    RV_dist = np.array([np.array(RV_dist[0]).astype(float),np.array(RV_dist[1]).astype(float)])
    temp_file = open('RV_dist_histogram'+args.plot_suffix+'.pkl', 'w')
    pickle.dump((RV_dist), temp_file)
    temp_file.close()

plt.clf()
no_bins = 30
region_labels = ['Upper Scorpius', 'Upper Centaurus-Lupus']
colors = ['b', 'orange']
sp_label = ['F', 'G', 'K', 'M']
bins = np.logspace(-1.0, 2.0, no_bins)
centers = (bins[:-1] + bins[1:])/2.
hist_US, bins = np.histogram(RV_dist[0], bins=bins)
hist_UCL, bins = np.histogram(RV_dist[1][:57], bins=bins)
y = hist_US + hist_UCL
popt, pcov = optimize.curve_fit(function,centers, y)
x_fit = np.logspace(-1.0, 2.0, 1000)
y_fit = function(x_fit, *popt)

width = bins[1:] - bins[:-1]
h_US = plt.bar(centers, hist_US, width=width, align='center', label='Upper Scorpius', color='b')
h_UCL = plt.bar(centers, hist_UCL, width=width, align='center', label='Upper Centaurus-Lupus', bottom=hist_US, color='orange')
#plt.legend(loc='best')
plt.plot(x_fit, y_fit, c='k')
plt.title('mean:' + str(x_fit[np.argmax(y_fit)]) + ' std:' + str(popt[2]))
plt.xscale('log')
plt.xlabel('$\Delta RV$ (km/s)')
plt.ylim([0, 25.])
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Observational_Paper/Images/RV_distributions/Radial_Velocity_distribution' + args.plot_suffix + '_' + str(no_bins) + '.eps')

'''
for dist in range(len(RV_dist)):
    if np.max(np.max(RV_dist[dist])) > max:
        max = np.max(np.max(RV_dist[dist]))
    bins = np.arange(0, max+2,2)
    hist, bins = np.histogram(RV_dist[dist], bins=bins)
    width = 0.5 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2. - width/2.
    plt.bar(center+(width*dist), hist, align='center', width=width, label=region_labels[dist], color=colors[dist], alpha=0.5)

    pl.figure()
    pl.hist(np.array(RV_dist[0]).astype(float), bins=np.logspace(0.1, 2.0, 50),alpha=0.50)
    pl.hist(np.array(RV_dist[1]).astype(float), bins=np.logspace(0.1, 2.0, 50),alpha=0.50)
    pl.gca().set_xscale("log")
    pl.show()
'''

plt.legend(loc='best')
plt.xlabel('RV Variation (km/s)')
plt.ylabel('#')
plt.xlim([0.0, 100.0])
plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/Radial_Velocity_distribution' + args.plot_suffix + '.pdf')
print "CREATED RV DISTRIBUTION HISTOGRAM"

for Spt in range(len(RV_dist_sptype)):
    if len(RV_dist_sptype[Spt]) > 0:
        max = np.max(RV_dist_sptype[Spt])
        bins = np.arange(0, max+2,2)
        hist, bins = np.histogram(RV_dist_sptype[Spt], bins=bins)
        width = (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:])/2. - width/2.
        plt.clf()
        plt.bar(center+(width*dist), hist, align='center', width=width)
        pl.gca().set_xscale("log")
        plt.xlabel('RV Variation (km/s)')
        plt.ylabel('#')
        plt.xlim([0.0, 100.0])
        plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/Radial_Velocity_distribution_' + sp_label[Spt] + args.plot_suffix +'.pdf')
        print "CREATED RV DISTRIBUTION HISTOGRAM FOR SPECTRAL TYPE", sp_label[Spt]
    
    
    else:
        print "Couldn't make RV DISTRIBUTION HISTOGRAM FOR SPECTRAL TYPE", sp_label[Spt]


images_RIK = images = glob.glob('/Users/rajikak/Observational_Data/PDF_dirs/RIK*/RV_vs_MJD.png')
images_RIK.sort(key=lambda f: int(filter(str.isdigit, f)))
images_UCAC = images = glob.glob('/Users/rajikak/Observational_Data/PDF_dirs/UCAC*/RV_vs_MJD.png')
images_UCAC.sort(key=lambda f: int(filter(str.isdigit, f)))
images = images_RIK + images_UCAC
pdf = FPDF()
pdf.add_page()

x,y,w,h,row,column = 10,10,100,68,0,0
trim = 0.2
for image in images:
    if '.pdf' not in image:
        if (y+h) > 300:
            pdf.add_page()
            x,y,row,column = 10,10,0,0
        
        pdf.image(image, x, y, w, h)
        column = column + 1
        x = x + w
        if column == 2:
            column = 0
            x = 10
            y = y + h

pdf.output('/Users/rajikak/Observational_Data/PDF_dirs/Radial_Velocity_vs_MJD_all.pdf', 'F')
pdf.close()


if args.mode == 'make_paper_plot':
    #Histogram plot
    hist_pickles = ['RV_dist_histogram_highest_none.pkl', 'RV_dist_histogram_highest_abs_all.pkl']
    text_labels = ['No correction to retrieved RV', 'Correction using', 'B-band absorption']
    
    plt.clf()
    no_bins = 30
    region_labels = ['Upper Scorpius', 'Upper Centaurus-Lupus']
    colors = ['b', 'orange']
    
    bins = np.logspace(-1.0, 2.0, no_bins)
    centers = (bins[:-1] + bins[1:])/2.
    width = bins[1:] - bins[:-1]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    fig.set_size_inches(7,9)
    text_ypos = 23
    
    RV_dist = pickle.load(open(hist_pickles[0], "rb" ))
    hist_US, bins = np.histogram(RV_dist[0], bins=bins)
    hist_UCL, bins = np.histogram(RV_dist[1][:57], bins=bins)
    y = hist_US + hist_UCL
    popt, pcov = optimize.curve_fit(function,centers, y)
    x_fit = np.logspace(-1.0, 2.0, 1000)
    y_fit = function(x_fit, *popt)
    

    h_US = ax1.bar(centers, hist_US, width=width, align='center', label='Upper Scorpius', color='b')
    h_UCL = ax1.bar(centers, hist_UCL, width=width, align='center', label='Upper Centaurus-Lupus', bottom=hist_US, color='orange')
    ax1.plot(x_fit, y_fit, c='k')
    ax1.axvline(x = x_fit[np.argmax(y_fit)], linestyle='--', color='k')
    ax1.set_xscale('log')
    ax1.set_ylim([0, 25.])
    ax1.set_ylabel('#')
    ax1.text(1.5e-1, text_ypos, text_labels[0])
    print "median RV variation =", np.median(np.concatenate((RV_dist[0], RV_dist[1])))

    RV_dist = pickle.load(open(hist_pickles[1], "rb" ))
    hist_US, bins = np.histogram(RV_dist[0], bins=bins)
    hist_UCL, bins = np.histogram(RV_dist[1][:57], bins=bins)
    y = hist_US + hist_UCL
    popt, pcov = optimize.curve_fit(function,centers, y)
    x_fit = np.logspace(-1.0, 2.0, 1000)
    y_fit = function(x_fit, *popt)
    
    
    h_US = ax2.bar(centers, hist_US, width=width, align='center', label='Upper Scorpius', color='b')
    h_UCL = ax2.bar(centers, hist_UCL, width=width, align='center', label='Upper Centaurus-Lupus', bottom=hist_US, color='orange')
    ax2.plot(x_fit, y_fit, c='k')
    ax2.axvline(x = x_fit[np.argmax(y_fit)], linestyle='--', color='k')
    ax2.set_xscale('log')
    ax2.set_ylim([0, 25.])
    ax2.set_xlabel('$\Delta v_r$ (km/s)')
    ax2.set_ylabel('#')
    ax2.legend(loc='upper right')
    ax2.text(1.5e-1, text_ypos, text_labels[1])
    ax2.text(1.5e-1, text_ypos-1.5, text_labels[2])
    print "median RV variation =", np.median(np.concatenate((RV_dist[0], RV_dist[1])))
    plt.setp([ax2.get_yticklabels()[-1]], visible=False)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.subplots_adjust(hspace=0)
    plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Observational_Paper/Images/RV_distributions/Comp_Radial_Velocity_distribution' + args.plot_suffix + '_' + str(no_bins) + '.eps', format='eps', bbox_inches='tight')
    
    plt.clf()
    fig.set_size_inches(7,5)
    h_US = plt.bar(centers, hist_US, width=width, align='center', label='Upper Scorpius', color='b')
    h_UCL = plt.bar(centers, hist_UCL, width=width, align='center', label='Upper Centaurus-Lupus', bottom=hist_US, color='orange')
    plt.plot(x_fit, y_fit, c='k')
    plt.xscale('log')
    plt.ylim([0, 25.])
    plt.xlabel('$\Delta v_{r,extrema}$ (km/s)')
    plt.ylabel('#')
    plt.legend(loc='upper right')
    #plt.text(1.5e-1, text_ypos, text_labels[1])
    #plt.text(1.5e-1, text_ypos-1.5, text_labels[2])
    print "median RV variation =", np.median(np.concatenate((RV_dist[0], RV_dist[1])))
    plt.setp([ax2.get_yticklabels()[-1]], visible=False)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.subplots_adjust(hspace=0)
    plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Observational_Paper/Images/RV_distributions/Radial_Velocity_distribution' + args.plot_suffix + '_' + str(no_bins) + '.eps', format='eps', bbox_inches='tight')
    

    #Cross-correlation plots
    RIK_96_pickles = ['/Users/rajikak/Observational_Data/PDF_dirs/RIK-96/cross_correlation_56824.6915456.pkl', '/Users/rajikak/Observational_Data/PDF_dirs/RIK-96/cross_correlation_57498.6251913.pkl']
    
    plt.clf()
    fig.clf()
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    fig.set_size_inches(11,5)
    x,y = pickle.load(open(RIK_96_pickles[0], "rb" ))
    ax1.plot(x, y)
    #ax1.set_ylim(max=1.)
    ax1.set_xlim([-800, 800])
    ax1.set_ylabel('Cross-correlation fit')
    ax1.set_xlabel('Velocity (km/s)')
    ax1.text(-725, 1.05, 'RIK-96')

    x,y = pickle.load(open(RIK_96_pickles[1], "rb" ))
    ax2.plot(x, y/np.max(y))
    #ax2.set_ylim(max=1.)
    ax2.set_xlim([-800, 800])
    ax2.set_xlabel('Velocity (km/s)')
    plt.setp([ax2.get_yticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([ax1.get_xticklabels()[-1]], visible=False)
    fig.subplots_adjust(wspace=0)
    plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Observational_Paper/Images/SB2/RIK-96.eps', format='eps', bbox_inches='tight')

    #Cross-correlation plots
    UCAC4_pickles = ['/Users/rajikak/Observational_Data/PDF_dirs/UCAC4-161328427/cross_correlation_57477.5794273.pkl', '//Users/rajikak/Observational_Data/PDF_dirs/UCAC4-161328427/cross_correlation_57908.508085.pkl']
    
    plt.clf()
    fig.clf()
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    fig.set_size_inches(11,5)
    x,y = pickle.load(open(UCAC4_pickles[0], "rb" ))
    ax1.plot(x, y)
    #ax1.set_ylim(max=1.)
    ax1.set_xlim([-800, 800])
    ax1.set_ylabel('Cross-correlation fit')
    ax1.set_xlabel('Velocity (km/s)')
    ax1.text(-725, 1.1, 'UCAC4-161328427')
    
    x,y = pickle.load(open(UCAC4_pickles[1], "rb" ))
    ax2.plot(x, y/np.max(y))
    #ax2.set_ylim(max=1.)
    ax2.set_xlim([-800, 800])
    ax2.set_xlabel('Velocity (km/s)')
    plt.setp([ax2.get_yticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([ax1.get_xticklabels()[-1]], visible=False)
    fig.subplots_adjust(wspace=0)
    plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Observational_Paper/Images/SB2/UCAC4-161328427.eps', format='eps', bbox_inches='tight')