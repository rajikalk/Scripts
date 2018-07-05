import process_stellar as ps
import glob
import astropy.io.fits as pyfits
import csv
import numpy as np
import pickle
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
from scipy import optimize
import pylab as py
import os

def function(x, a, b, c):
    return a/(np.sqrt(2*np.pi)*c) * py.exp(-(x - b)**2.0 / (2 * c**2))

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-abs", "--absorption", help="do you want to correct for sky absorption?", default='False')
    parser.add_argument("-allspt", "--all_spectral_types", help="do you want to correct absorption for all spectral types", default='False')
    parser.add_argument("-sky", "--skylines", help="do you want to correct for skylines?", default='False')
    parser.add_argument("-abs_spt", "--abs_corr_spectral_types", help="Which spectral types would you like to correct the absorption line rv for? default is 'F,G'", type=str, default='F,G')
    parser.add_argument("-temp_dir", "--template_dir", help="what template directory do you want to use? If none, it uses the default in the /tools directory", type=str, default=None)
    parser.add_argument("-plot_suffix", "--plot_suffix", help="Do you want to give your plots suffixes to help distinguish between them?", type=str, default='')
    parser.add_argument("-mpp", "--make_paper_plot", type=str, default='False')
    parser.add_argument("files", nargs='*')
    
    args = parser.parse_args()
    return args

args = parse_inputs()

pickle_files = ['histogram_correction_none.pkl', 'histogram_correction_abs.pkl', 'histogram_correction_abs_FG.pkl']

reduced_dir = '/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/'

#ADD IN SECTION REMAKING THE TEMPLATES FOR EACH CONFIGURATION OF WHETHER TO USE SKY LINES OR ABSORPTION

#Read in Standards Data
RV_standard_name = []
RV_standard_data = []
RV_obj = []
RV_data = []
RV_abs = []
RV_abs_sig = []
Hcorr = []
D_rvs = [[],[],[],[]]
if args.template_dir == None:
    templates_dir = '/Users/rajikak/tools/templates/'
else:
    templates_dir = args.template_dir
abs_temps = glob.glob('/Users/rajikak/tools/absorption_spec_F.fits')
sky_temp = glob.glob('/Users/rajikak/tools/wifes_sky.fits')
abs_corr_spt = args.abs_corr_spectral_types.split(',')

sky_intervals = ([0,5500],[5700,5850],[6000,6100],[6700,6800])
abs_bad_ints = ([0,5500],[5500,6860],[6910,7000])

header = 0
with open('/Users/rajikak/Observational_Data/RV_standards/RV_standard_clean.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            RV_standard_name.append(row[0])
            RV_standard_data.append(row[1:])
        else:
            header = 1

f = open('/Users/rajikak/Observational_Data/RV_standards/RV_standard_testing.csv', 'w')
f.write('Name,RV(km/s),Spectral Type,RA,DEC,file,SNR,Retrieved RV,RV_sig,Delta_RV,template\n')

if args.make_paper_plot == 'False':
    images = glob.glob('/Users/rajikak/Observational_Data/RV_standards/Images/*/*fits*')
    for image in images:
        obj_name = image.split('/')[-2]
        if obj_name in RV_standard_name:
            ind = RV_standard_name.index(obj_name)
            file_dir = image.split('/')[-1].split('-')[1].split('.')[0]
            file = image.split('/')[-1].split('.png')[0]

            full_file_path = reduced_dir + file_dir + '/' + file
            
            if os.path.isfile(full_file_path):
                hdu = pyfits.open(full_file_path)
                h_corr = hdu[0].header['RADVEL']
                Hcorr.append(h_corr)

                flux,wave,var,sky = ps.read_and_find_star_p08(full_file_path, sky_rad=5, return_sky=True)
                spectrum,sig = ps.weighted_extract_spectrum(flux,var)
                
                if args.absorption != 'False':
                    if args.all_spectral_types != 'False':
                        rv_abs, rv_abs_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,abs_temps, (abs_bad_ints), heliocentric_correction=0.0)
                        print "ABSORPTION OFFSET= " + str(rv_abs) + ', '+ str(rv_abs_sig)
                    else:
                        if RV_standard_data[ind][1][0] in abs_corr_spt:
                            rv_abs, rv_abs_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,abs_temps, (abs_bad_ints), heliocentric_correction=0.0)
                            print "ABSORPTION OFFSET= " + str(rv_abs) + ', '+ str(rv_abs_sig)
                        else:
                            rv_abs = 0.0
                            rv_abs_sig = np.nan
                    if np.isnan(rv_abs) or np.isnan(rv_abs_sig):
                        print "COULD NOT GET A RELIABLE ABSORPTION RV"
                        rv_abs = 0.0
                    if np.isnan(rv_abs_sig):
                        rv_abs_sig = np.inf
                else:
                    rv_abs = 0.0

                if args.skylines != 'False':
                    spectrum,sig = ps.weighted_extract_spectrum(sky,var)
                    rv_sky, rv_sky_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,sky_temp, (sky_intervals), heliocentric_correction=0.0)
                    print "SKYLINE OFFSET= " + str(rv_sky) + ', '+ str(rv_sky_sig)
                    if np.isnan(rv_sky) or np.isnan(rv_sky_sig) or np.abs(rv_sky)>100.:
                        print "COULD NOT GET A RELIABLE SKYLINE RV"
                        rv_sky = 0.0
                        rv_sky_sig = np.inf
                    if np.isnan(rv_sky_sig):
                        rv_sky_sig = np.inf
                else:
                    rv_sky = 0.0

                if args.absorption != 'False' and args.skylines != 'False':
                    if rv_sky_sig < rv_abs_sig:
                        rv_abs = 0.0
                        print("USING CORRECTION FROM SKYLINES")
                    elif rv_sky_sig > rv_abs_sig:
                        rv_sky = 0.0
                        print("USING CORRECTION FROM ABSORPTION")


                template_name = '/Users/rajikak/tools/templates/' + obj_name + '_' + file.split('.p08')[0]+'.fits'

                usable_temps = glob.glob(templates_dir + '*T2m3w*.fits')
                if templates_dir + template_name.split('/')[-1] in usable_temps:
                    usable_temps.remove(templates_dir + template_name.split('/')[-1])

                spectrum,sig = ps.weighted_extract_spectrum(flux,var)
                rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,usable_temps, ([0,5400],[6900,7000]), save_figures=False, heliocentric_correction=(h_corr - rv_abs -rv_sky))
                RV_data.append(rv)

                #rv = rv - rv_abs

                custom_simbad = Simbad()
                custom_simbad.add_votable_fields("sptype")
                data = custom_simbad.query_object(obj_name)

                Delta_RV = float(RV_standard_data[ind][0]) - rv

                sptype = data['SP_TYPE'][-1]
                if (sptype[0] == 'F'):
                    D_rvs[0].append(Delta_RV)
                if (sptype[0] == 'G'):
                    D_rvs[1].append(Delta_RV)
                if (sptype[0] == 'K'):
                    D_rvs[2].append(Delta_RV)
                if (sptype[0] == 'M'):
                    D_rvs[3].append(Delta_RV)

                print "######################################"
                print "#"
                print "#"
                print "#"
                print "#"
                print "FILE =", obj_name + '_' + file
                print "TEMP_USED", temp_used
                print "SPECTRAL TYPE =", sptype
                print "OBJECT RV =", RV_standard_data[ind][0]
                print "OBTAINED RV=", rv
                print "DELTA RV=", np.abs(float(RV_standard_data[ind][0]) - float(rv))
                print "#"
                print "#"
                print "#"
                print "#"
                print "######################################"


                write_line = obj_name + ',' + ','.join(RV_standard_data[ind]) + ',' + file + ',' + str(np.mean(spectrum)) + ',' + str(rv) + ',' + str(rv_sig) + ',' + str(Delta_RV) + ',' + temp_used + '\n'
                f.write(write_line)
            else:
                print "SKIPPING IMAGE"
                continue

    f.close()

    bins = np.arange(np.floor(np.hstack(D_rvs).min()), np.ceil(np.hstack(D_rvs).max())+1)
    #hist_F, bins = np.histogram(D_rvs[0], bins=bins)
    hist_G, bins = np.histogram(D_rvs[1], bins=bins)
    hist_K, bins = np.histogram(D_rvs[2], bins=bins)
    hist_M, bins = np.histogram(D_rvs[3], bins=bins)
    width = 1
    centers = (bins[:-1] + bins[1:])/2.
    y = hist_G + hist_K + hist_M
    popt, pcov = optimize.curve_fit(function,centers, y)
    x_fit = np.linspace(bins[0], bins[-1], 100)
    y_fit = function(x_fit, *popt)
    total_gauss_fit = (x_fit, y_fit)

    plt.clf()
    h_M = plt.bar(centers, hist_M, align='center', width=width, label='M', color='m')
    h_K = plt.bar(centers, hist_K, align='center', width=width, label='K', color='r', bottom=hist_M)
    h_G = plt.bar(centers, hist_G, align='center', width=width, label='G', color='g', bottom=(hist_M+hist_K))
    #h_F = plt.bar(centers, hist_F, align='center', width=width, label='F', color='b', bottom=(hist_M+hist_K+hist_G))
    plt.plot(x_fit, y_fit)
    handles = [h_M, h_K, h_G]#, h_F]
    plt.legend(loc='best', handles=handles[::-1])
    plt.xlabel('$\Delta v_r$ (km/s)')
    plt.ylabel('#')
    plt.xlim([-8,8])
    plt.ylim([0,20])
    #plt.title('RV precision (mu:' + np.str(np.mean(np.hstack(D_rvs))) +', sigma:' + np.str(np.std(np.hstack(D_rvs), ddof=1)) + ')')
    plt.savefig('RV_precision'+args.plot_suffix+'.eps')

    '''
    y = hist_F
    popt_F, pcov_F = optimize.curve_fit(function,centers, y)
    x_fit = np.linspace(bins[0], bins[-1], 100)
    y_fit = function(x_fit, *popt_F)
    '''

    y = hist_G
    popt_G, pcov_G = optimize.curve_fit(function,centers, y)
    x_fit = np.linspace(bins[0], bins[-1], 100)
    y_fit = function(x_fit, *popt_G)

    y = hist_K
    popt_K, pcov_K = optimize.curve_fit(function,centers, y)
    x_fit = np.linspace(bins[0], bins[-1], 100)
    y_fit = function(x_fit, *popt_K)

    y = hist_M
    #popt_M, pcov_M = optimize.curve_fit(function,centers, y)
    #x_fit = np.linspace(bins[0], bins[-1], 100)
    #y_fit = function(x_fit, *popt_M)

    #print "For F: mu = ", np.mean(D_rvs[0]), "std = ", np.std(D_rvs[0], ddof=1)
    print "For G: mu = ", np.mean(D_rvs[1]), "std = ", np.std(D_rvs[1], ddof=1)
    print "For K: mu = ", np.mean(D_rvs[2]), "std = ", np.std(D_rvs[2], ddof=1)
    print "For M: mu = ", np.mean(D_rvs[3]), "std = ", np.std(D_rvs[3], ddof=1)
    print "Total: mu = ", popt[1], "std = ", popt[2]


    temp_file = open('histogram'+args.plot_suffix+'.pkl', 'w')
    pickle.dump((centers, hist_M, hist_K, hist_G, total_gauss_fit), temp_file)
    temp_file.close()


if args.make_paper_plot == 'True':
    plt.clf()
    fig = plt.figure(figsize=(6,10))
    text_labels = ['Interpolation between arcs', 'for all standards', 'for G standards']
    for plot_pickle in range(len(pickle_files)):
        centers, hist_M, hist_K, hist_G, total_gauss_fit = pickle.load(open(pickle_files[plot_pickle], "rb" ))
        text_ypos = 12.5
        if plot_pickle == 0:
            ax1 = fig.add_subplot(3,1,plot_pickle+1)
            ax1.set_xlim([-8,8])
            ax1.set_ylim([0,15])
            h_M = ax1.bar(centers, hist_M, align='center', width=1, label='M', color='r')
            h_K = ax1.bar(centers, hist_K, align='center', width=1, label='K', color='g', bottom=hist_M)
            h_G = ax1.bar(centers, hist_G, align='center', width=1, label='G', color='b', bottom=(hist_M+hist_K))
            #h_F = ax1.bar(centers, hist_F, align='center', width=1, label='F', color='b', bottom=(hist_M+hist_K+hist_G))
            ax1.plot(total_gauss_fit[0], total_gauss_fit[1])
            handles = [h_M, h_K, h_G]#, h_F]
            ax1.legend(loc='best', handles=handles[::-1])
            ax1.set_ylabel('#')
            ax1.text(-7.5, text_ypos, text_labels[plot_pickle])
        else:
            ax = fig.add_subplot(3,1,plot_pickle+1, sharex=ax1, sharey=ax1)
            h_M = ax.bar(centers, hist_M, align='center', width=1, label='M', color='r')
            h_K = ax.bar(centers, hist_K, align='center', width=1, label='K', color='g', bottom=hist_M)
            h_G = ax.bar(centers, hist_G, align='center', width=1, label='G', color='b', bottom=(hist_M+hist_K))
            #h_F = ax.bar(centers, hist_F, align='center', width=1, label='F', color='b', bottom=(hist_M+hist_K+hist_G))
            ax.plot(total_gauss_fit[0], total_gauss_fit[1])
            ax.set_ylabel('#')
            text_ypos = 12.5
            plt.setp([ax.get_yticklabels()[-1]], visible=False)
            ax.text(-7.5, text_ypos, 'B-band correction')
            ax.text(-7.5, text_ypos-1.5, text_labels[plot_pickle])
    ax.set_xlabel('$v_{r,standard} -  v_{r,retrieved}$ (km/s)')
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    fig.subplots_adjust(hspace=0)
    fig.savefig('RV_precision.eps', format='eps', bbox_inches='tight')