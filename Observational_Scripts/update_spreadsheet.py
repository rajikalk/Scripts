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
RV_variation = []
Temp_sptype = []
Pref_template = []
Obs_info = []
SNR = []

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
            RV_variation.append(row[13])
            Pref_template.append(row[14])
            Temp_sptype.append(row[15])
            if len(row) > 16:
                Obs = np.array(row[16:])
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

#Reduced dirs
if args.mode == 'update':
    p08_dirs = glob.glob('Observational_Data/Reduced_Wifes_Obs/2*/')
    if args.start_d != '':
        start_ind = p08_dirs.index('Observational_Data/Reduced_Wifes_Obs/'+args.start_d+'/')
    else:
        start_ind = 0
    if args.end_d != '':
        end_ind = p08_dirs.index('Observational_Data/Reduced_Wifes_Obs/'+args.end_d+'/')
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
            Obs_date = hdu[0].header['DATE-OBS'].split('T')[0]
            MJD = hdu[0].header['MJD-OBS']
            #alt_date = Obs_date.split('-')[-1]+'/'+Obs_date.split('-')[-2]+'/'+Obs_date.split('-')[-3][2:]
            
            if [Obj_name, Obs_date] in ignore_obs_list:
                print "SKIPPING FILE:", fn, " BECAUSE DATA IS NOT GOOD"
            elif (Obj_name in Object and args.obj_name == "") or (Obj_name == args.obj_name) or [item for item in Manual_click_flag if item[0] == Obj_name and item[1] == Obs_date] != []:
                ind = Object.index(Obj_name)
                if (str(np.round(MJD, decimals=5)) not in Obs_info[ind]) or (Obs_date == args.date) or ((Obj_name, Obs_date) in Manual_click_flag):
                    h_corr = hdu[0].header['RADVEL']
                    RA = hdu[0].header['RA']
                    DEC = hdu[0].header['DEC']
                    
                    if args.click == 'True' or [item for item in Manual_click_flag if item[0] == Obj_name and item[1] == Obs_date] != []:
                        click = True
                        manual_click_data = [item for item in Manual_click_flag if item[0] == Obj_name and item[1] == Obs_date]
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
                       
                       if np.isnan(rv_abs) or np.isnan(rv_abs_sig): # or Temp_sptype[ind][0] == 'M' or Temp_sptype[ind][0] == 'K':
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
                        spectrum,sig = ps.weighted_extract_spectrum(flux)
                        rv_sky = 0.0
                        rv_abs = 0.0
                        rv_sky_sig = 0.0
                        rv_abs_sig = 0.0
                        rv_err = 0.0
                        which_sky = ''
                    
                    #selecting template
                    if Pref_template[ind] != '':
                        templates = glob.glob('/Users/rajikak/tools/templates/'+Pref_template[ind]+'*')
                        if Temp_sptype[ind] == '':
                            temp_name = Pref_template[ind].split('_')[0]
                            custom_simbad = Simbad()
                            custom_simbad.add_votable_fields("sptype")
                            data = custom_simbad.query_object(temp_name)
                            sptype = data['SP_TYPE'][-1]
                            Temp_sptype[ind] = sptype
                        else:
                            sptype = Temp_sptype[ind]
                        '''
                        if sptype[0] == 'M':
                            rv_abs = 0.0
                            rv_sky = 0.0
                        else:
                            which_sky = 'abs'
                        '''
                    else:
                        templates = glob.glob('/Users/rajikak/tools/templates/*.fits')

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
RV_dist = [[[],[]],[[],[]]]
RV_dist_sptype = [[], [], [], []]
#RV_dist = []
for obj in range(len(Object)):
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
            if SB2_flag[obj] == True or SB1_flag[obj] == True:
                RV_dist[0][1].append(spread_val)
                print "2: Added rv spread of a binary US disk bearing object"
            else:
                RV_dist[0][0].append(spread_val)
                print "1: Added rv spread of a US disk bearing object"
        elif 'Y' in IR_excess[obj]:
            if SB2_flag[obj] == True or SB1_flag[obj] == True:
                RV_dist[1][1].append(spread_val)
                print "4: Added rv spread of a binary UCL disk bearing object"
            else:
                RV_dist[1][0].append(spread_val)
                print "3: Added rv spread of a UCL disk bearing object"

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
    images = glob.glob('Observational_Data/PDF_dirs/'+Object[obj]+'/*')
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

    pdf.output('Observational_Data/PDF_dirs/'+Object[obj]+'/'+Object[obj]+'_all.pdf', 'F')
    pdf.close()

#write out updated data
f = open(args.input_file, 'w')
f.write('Object,Region,RA,DEC,Pmem,Disk,Visual inspection,K_mag,V_mag,B_mag,No. Obs,SB1_flag,SB2 flag,RV_variation,Pref_temp,Temp_SpT,Observations,MJD,RV,RV_err,RV_sky,Template\n')

for obj in range(len(Object)):
    write_string = Object[obj] + ',' + Region[obj] + ',' + Coords[obj][0] + ',' + Coords[obj][1] + ',' + str(Membership_probability[obj]) + ',' + IR_excess[obj] + ',' + Disk_visual[obj] + ',' + str(Magnitudes[obj][0]) + ',' + str(Magnitudes[obj][1]) + ',' + str(Magnitudes[obj][2]) + ',' + str(No_obs[obj]) + ',' + str(SB1_flag[obj]) + ','+ str(SB2_flag[obj]) + ',' +str(RV_variation[obj]) + ',' + Pref_template[obj] + ',' + Temp_sptype[obj] + ',' + ','.join(np.ravel(Obs_info[obj])) + '\n'
    f.write(write_string)

f.close()

plt.clf()
max = 0
region_labels = ['Upper Scorpius', 'Upper Centaurus-Lupus']
colors = ['b', 'orange']
sb_colors = ['g', 'magenta']
sp_label = ['F', 'G', 'K', 'M']
for dist in range(len(RV_dist)):
    if np.max(np.max(RV_dist[dist])) > max:
        max = np.max(np.max(RV_dist[dist]))
    bins = np.arange(0, max+2,2)
    hist_s, bins_s = np.histogram(RV_dist[dist][0], bins=bins)
    hist_b, bins_b = np.histogram(RV_dist[dist][1], bins=bins)
    width = 0.5 * (bins_s[1] - bins_s[0])
    center = (bins_s[:-1] + bins_s[1:])/2. - width/2.
    plt.bar(center+(width*dist), hist_s, align='center', width=width, label=region_labels[dist], color=colors[dist], alpha=0.5)
    plt.bar(center+(width*dist), hist_b, align='center', width=width, label=region_labels[dist]+'(SB)', color=sb_colors[dist], hatch='//', bottom=hist_s)

plt.legend(loc='best')
plt.xlabel('RV Variation (km/s)')
plt.ylabel('#')
plt.xlim([0.0, 100.0])
plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/Radial_Velocity_distribution' + args.plot_suffix + '.pdf')
print "CREATED RV DISTRIBUTION HISTOGRAM"

for Spt in range(len(RV_dist_sptype)):
    max = np.max(RV_dist_sptype[Spt])
    bins = np.arange(0, max+2,2)
    hist, bins = np.histogram(RV_dist_sptype[Spt], bins=bins)
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2. - width/2.
    plt.clf()
    plt.bar(center+(width*dist), hist, align='center', width=width)
    plt.xlabel('RV Variation (km/s)')
    plt.ylabel('#')
    plt.xlim([0.0, 100.0])
    plt.savefig('/Users/rajikak/Observational_Data/PDF_dirs/Radial_Velocity_distribution_' + sp_label[Spt] + args.plot_suffix +'.pdf')
    print "CREATED RV DISTRIBUTION HISTOGRAM FOR SPECTRAL TYPE", sp_label[Spt]


images_RIK = images = glob.glob('Observational_Data/PDF_dirs/RIK*/RV_vs_MJD.png')
images_RIK.sort(key=lambda f: int(filter(str.isdigit, f)))
images_UCAC = images = glob.glob('Observational_Data/PDF_dirs/UCAC*/RV_vs_MJD.png')
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

pdf.output('Observational_Data/PDF_dirs/Radial_Velocity_vs_MJD_all.pdf', 'F')
pdf.close()


