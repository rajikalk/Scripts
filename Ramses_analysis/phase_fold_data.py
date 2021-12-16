import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os
import yt
from scipy.stats import mode
import numpy.ma as ma
from scipy.stats import skewnorm
import scipy.stats as stats
from scipy import optimize
from scipy.integrate import simps
from scipy.optimize import curve_fit
import math as math
import scipy.special as sp
import numpy.ma as ma

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-def_sys", "--define_system", help="Is there a particular system that you would like to plot?", type=str, default=None)
    parser.add_argument("-tf", "--text_font", help="text font for figure", default=12, type=int)
    parser.add_argument("-pd", "--pickle_dump", help="do you want to pickle the phasefolded data?", default='False', type=str)
    parser.add_argument("-e_bins", "--eccentricity_bins", help="define bins", default='[1.1,0.6,0.5,0.4,0.3,0.2,0.1,0.0]', type=str)
    parser.add_argument("-n_bin", "--number_of_bins", help="Default number of bins is 20", default=20, type=int)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]
    
#================================================================================
args = parse_inputs()

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)
    
systems_hierarchy = args.define_system

#read in particle data:
try:
    sytem_pickle = glob.glob(path + 'reduced_system_data.pkl')[0]
    file_open = open(sytem_pickle, 'rb')
    reduced_systems_data = pickle.load(file_open)
    file_open.close()
    time = reduced_systems_data['time']
    print("Finished reading in reduced systems data")
except:
    print("reduced_system_data.pkl doesn't exist. Genrate it first using binary_evolution_plot.py")
    
#Create loop and load binary data:
system_tag_int = 65
if systems_hierarchy != None:
    while len(find(systems_hierarchy, '[')) > 0:
        indexes = []
        for ind, char in enumerate(systems_hierarchy):
            print('systems_hierarchy =', systems_hierarchy)
            if char == '[':
                indexes.append(ind)
            if char == ']':
                start_ind = indexes.pop()
                end_ind = ind
                
                binary_tags = eval('['+systems_hierarchy[start_ind+1:end_ind]+']')
                system_tag = chr(system_tag_int)
                system_tag_int = system_tag_int + 1
                       
                systems_hierarchy = systems_hierarchy[:start_ind] + "\'" + system_tag + "\'" + systems_hierarchy[end_ind+1:]
                
                system_ind = reduced_systems_data['base_tags'].index(binary_tags)
                print("CHECK SEPARATION CALCULATION TAKES INTO ACCOUNT CYCLIC BOUNDARIES")
                separation = reduced_systems_data['separation'][system_ind]
                eccentricity = reduced_systems_data['eccentricity'][system_ind]
                accretion = reduced_systems_data['mdot_individual'][system_ind]
                
                ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
                ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
                periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
                apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]
                
                if len(periastron_inds) > 10:
                    e_bins = eval(args.eccentricity_bins)
                    e_bin_it = 1
                    plot_it = 0
                    phase = np.linspace(0.0, 1.0, args.number_of_bins+1)
                    phase_centers = (phase[1:] + phase[:-1])/2.
                    phase_pre = np.linspace(-1.0, 0.0, args.number_of_bins+1)
                    phase_2 = np.linspace(1.0, 2.0, args.number_of_bins+1)
                    phase_2 = phase_pre.tolist()[:-1] + phase.tolist() + phase_2[1:].tolist()
                    phase_centers = (np.array(phase_2[1:]) + np.array(phase_2[:-1]))/2.
                    
                    median_eccentricity = []
                    std_eccentricity = []
                    
                    multiple_folds = []
                    accretion_err = []
                    y_fits = []
                    
                    plt.clf()
                    fig = plt.figure()
                    total_number_of_figures = len(e_bins) - 1
                    columns = 5
                    rows = int(np.ceil(total_number_of_figures/columns))
                    fig.set_size_inches(3.25*columns, 3.25*rows)
                    gs = gridspec.GridSpec(rows, columns)
                    
                    gs.update(wspace=0.0, hspace=0.0)
                    axes_dict = {}
                    
                    while e_bin_it < len(e_bins):
                        file_name = save_dir + str(binary_tags[0]) + "_" + str(binary_tags[1]) + '_accretion_median_start_orbit_from_' + str(e_bins[e_bin_it-1]) + '_' + str(e_bins[e_bin_it])
                        usable_periastrons = np.argwhere((eccentricity[periastron_inds]<e_bins[e_bin_it-1])&(eccentricity[periastron_inds]>e_bins[e_bin_it])).T[0]
                        if len(usable_periastrons) <2:
                            e_bin_it = e_bin_it + 1
                        else:
                            median_e = np.median(eccentricity[periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
                            std_e = np.std(eccentricity[periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
                            mean_e = np.mean(eccentricity[periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
                            median_eccentricity.append(median_e)
                            std_eccentricity.append([median_e-(mean_e-std_e), (mean_e+std_e)-median_e])
                            
                            rit = 1
                            averaged_binned_accretion = [[],[]]
                            averaged_total_accretion = []
                            print("e=["+str(e_bins[e_bin_it-1])+','+str(e_bins[e_bin_it])+"] --> "+str(len(usable_periastrons)) + " orbits")
                            
                            while rit < len(usable_periastrons):
                                binned_accretion = [[],[]]
                                total_accretion = []
                                
                                time_bins_1 = sorted(np.linspace(time[periastron_inds[usable_periastrons[rit-1]]], time[apastron_inds[usable_periastrons[rit-1]]],(int(args.number_of_bins/2)+1)))
                                time_bins_2 = sorted(np.linspace(time[apastron_inds[usable_periastrons[rit-1]]], time[periastron_inds[usable_periastrons[rit]]],(int(args.number_of_bins/2)+1)))
                                bin_ind = 1
                                
                                while bin_ind < len(time_bins_1):
                                    time_bin_inds_1 = np.where((time > time_bins_1[bin_ind-1]) & (time < time_bins_1[bin_ind]))[0]
                                    if len(time_bin_inds_1) != 0:
                                        intergrated_values = np.trapz(accretion[:,time_bin_inds_1], time[time_bin_inds_1])/(time[time_bin_inds_1[-1]]-time[time_bin_inds_1[0]])
                                        binned_accretion[0].append(intergrated_values[0])
                                        binned_accretion[1].append(intergrated_values[1])
                                        total_accretion.append(np.sum(np.nan_to_num(intergrated_values)))
                                    else:
                                        binned_accretion[0].append(yt.YTQuantity(np.nan, 'msun/yr'))
                                        binned_accretion[1].append(yt.YTQuantity(np.nan, 'msun/yr'))
                                        total_accretion.append(yt.YTQuantity(np.nan, 'msun/yr'))
                                    bin_ind = bin_ind + 1
                                bin_ind = 1
                                while bin_ind < len(time_bins_2):
                                    time_bin_inds_2 = np.where((time > time_bins_2[bin_ind-1]) & (time < time_bins_2[bin_ind]))[0]
                                    if len(time_bin_inds_2) != 0:
                                        intergrated_values = np.trapz(accretion[:,time_bin_inds_2], time[time_bin_inds_2])/(time[time_bin_inds_2[-1]]-time[time_bin_inds_2[0]])
                                        binned_accretion[0].append(intergrated_values[0])
                                        binned_accretion[1].append(intergrated_values[1])
                                        total_accretion.append(np.sum(np.nan_to_num(intergrated_values)))
                                    else:
                                        binned_accretion[0].append(yt.YTQuantity(np.nan, 'msun/yr'))
                                        binned_accretion[1].append(yt.YTQuantity(np.nan, 'msun/yr'))
                                        total_accretion.append(yt.YTQuantity(np.nan, 'msun/yr'))
                                    bin_ind = bin_ind + 1
                                    
                                averaged_binned_accretion[0].append(binned_accretion[0])
                                averaged_binned_accretion[1].append(binned_accretion[1])
                                averaged_total_accretion.append(total_accretion)
                                rit = rit + 1
                                
                            ax_label = 'ax' + str(plot_it)
                            if plot_it == 0:
                                axes_dict.update({ax_label:fig.add_subplot(gs[0,plot_it])})
                                axes_dict[ax_label].tick_params(axis="x",direction="in")
                                axes_dict[ax_label].set_xlim([0.0, 1.3])
                                #axes_dict[ax_label].set_ylim(bottom=0.0)
                                #axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                                axes_dict[ax_label].set_ylabel("Accretion Rate (M$_\odot$/yr)", fontsize=args.text_font)
                                axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
                                axes_dict[ax_label].tick_params(axis="y",direction="in")
                            else:
                                axes_dict.update({ax_label:fig.add_subplot(gs[int(plot_it/columns),np.remainder(plot_it,columns)], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
                                if plot_it < 1:
                                    xticklabels = axes_dict[ax_label].get_xticklabels()
                                    plt.setp(xticklabels, visible=False)
                                elif plot_it > 0:
                                    if plot_it != len(e_bins)-2:
                                        xticklabels = axes_dict[ax_label].get_xticklabels()
                                        plt.setp(xticklabels[0], visible=False)
                                elif plot_it > columns-1:
                                    if plot_it != len(e_bins)-2:
                                        xticklabels = axes_dict[ax_label].get_xticklabels()
                                        plt.setp(xticklabels[-2], visible=False)
                                if np.remainder(plot_it,columns) != 0:
                                    yticklabels = axes_dict[ax_label].get_yticklabels()
                                    plt.setp(yticklabels, visible=False)
                                    yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
                                    plt.setp(yticklabels, visible=False)
                                else:
                                    #axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                                    axes_dict[ax_label].set_ylabel("Accretion Rate (M$_\odot$/yr)", fontsize=args.text_font)
                                axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
                                axes_dict[ax_label].tick_params(axis="x",direction="in")
                                

                            median_accretion = []
                            median_accretion_total = []
                            standard_deviation = []
                            standard_deviation_total = []
                            
                            median_acc = np.median(averaged_binned_accretion[0], axis=0)
                            skip = True
                            try:
                                if True not in np.isnan(median_acc):
                                    skip = False
                            except:
                                if np.isnan(median_acc):
                                    skip = False
                            if skip == False:
                                mean_acc = np.mean(averaged_binned_accretion[0], axis=0)
                                std_acc = np.std(averaged_binned_accretion[0], axis=0)
                                median_accretion.append(median_acc)
                                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                                median_acc = np.median(averaged_binned_accretion[1], axis=0)
                                mean_acc = np.mean(averaged_binned_accretion[1], axis=0)
                                std_acc = np.std(averaged_binned_accretion[1], axis=0)
                                median_accretion.append(median_acc)
                                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                                
                                median_acc_tot = np.median(averaged_total_accretion, axis=0)
                                mean_acc_tot = np.mean(averaged_total_accretion, axis=0)
                                std_acc_tot = np.std(averaged_total_accretion, axis=0)
                                median_accretion_total.append(median_acc_tot)
                                standard_deviation_total.append([median_acc_tot-(mean_acc_tot-std_acc_tot), (mean_acc_tot+std_acc_tot)-median_acc_tot])
                                
                                long_median_accretion_1 = median_accretion[0].tolist()+median_accretion[0].tolist()+median_accretion[0].tolist()
                                long_median_accretion_2 = median_accretion[1].tolist()+median_accretion[1].tolist()+median_accretion[1].tolist()
                                long_median_total = median_accretion_total[0].tolist() + median_accretion_total[0].tolist() + median_accretion_total[0].tolist()
                                yerr_1 = [standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist(), standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()]
                                yerr_2 = [standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist(), standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()]
                                yerr = [standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist(), standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()]
                                accretion_err.append(np.array(yerr)/np.array(long_median_total))
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_1)*(1.e4), yerr=np.array(yerr_1)*(1.e4), ls='steps-mid', alpha=0.5, label='Primary')
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_2)*(1.e4), yerr=np.array(yerr_2)*(1.e4), ls='steps-mid', alpha=0.5, label='Secondary')
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_total)*(1.e4), yerr=np.array(yerr)*(1.e4), ls='steps-mid', label='Total')
                            else:
                                median_acc = np.median(np.nan_to_num(averaged_binned_accretion[0]), axis=0)
                                mean_acc = np.mean(np.nan_to_num(averaged_binned_accretion[0]), axis=0)
                                std_acc = np.std(np.nan_to_num(averaged_binned_accretion[0]), axis=0)
                                median_accretion.append(median_acc)
                                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                                median_acc = np.median(np.nan_to_num(averaged_binned_accretion[1]), axis=0)
                                mean_acc = np.mean(np.nan_to_num(averaged_binned_accretion[1]), axis=0)
                                std_acc = np.std(np.nan_to_num(averaged_binned_accretion[1]), axis=0)
                                median_accretion.append(median_acc)
                                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                                
                                median_acc_tot = np.median(averaged_total_accretion, axis=0)
                                mean_acc_tot = np.mean(averaged_total_accretion, axis=0)
                                std_acc_tot = np.std(averaged_total_accretion, axis=0)
                                median_accretion_total.append(mean_acc_tot)
                                standard_deviation_total.append([median_acc_tot-(mean_acc_tot-std_acc_tot), (mean_acc_tot+std_acc_tot)-median_acc_tot])
                                
                                long_median_accretion_1 = median_accretion[0].tolist()+median_accretion[0].tolist()+median_accretion[0].tolist()
                                long_median_accretion_2 = median_accretion[1].tolist()+median_accretion[1].tolist()+median_accretion[1].tolist()
                                long_median_total = median_accretion_total[0].tolist() + median_accretion_total[0].tolist() + median_accretion_total[0].tolist()
                                yerr_1 = [standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist(), standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()]
                                yerr_2 = [standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist(), standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()]
                                yerr = [standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist(), standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()]
                                accretion_err.append(np.array(yerr)/np.array(long_median_total))
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_1)*(1.e4), yerr=np.array(yerr_1)*(1.e4), ls='steps-mid', alpha=0.5, label='Primary')
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_2)*(1.e4), yerr=np.array(yerr_2)*(1.e4), ls='steps-mid', alpha=0.5, label='Secondary')
                                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_total)*(1.e4), yerr=np.array(yerr)*(1.e4), ls='steps-mid', label='Total')
                            time_text = plt.text(0.05, axes_dict[ax_label].get_ylim()[1]-(axes_dict[ax_label].get_ylim()[1] - axes_dict[ax_label].get_ylim()[0])*0.1, '$e$=['+str(e_bins[e_bin_it-1])+','+str(e_bins[e_bin_it])+'] using '+str(len(usable_periastrons))+' orbits', va="center", ha="left", color='k', fontsize=args.text_font)

                            #multiple_folds.append(np.array(long_median_total)*(1.e4))
                            multiple_folds.append(np.array(long_median_total))
                            '''
                            x_data = phase_centers[23:-15]
                            x = np.linspace(np.min(x_data),np.max(x_data),100)
                            if os.path.exists(file_name+'.pkl'):
                                file_open = open(file_name+'.pkl', 'rb')
                                phase_centers, long_median_accretion, yerr_tot, popt = pickle.load(file_open)
                                file_open.close()
                                y_fit= func(x, *popt)
                            else:
                                try:
                                    y_data = np.array(long_median_total)[23:-15]
                                    results = []
                                    pulsed_likelihoods = []
                                    no_pulsed_likelihoods = []
                                    for tries in range(1000):
                                        sigma = np.random.random()*2*0.2
                                        amp = np.random.random()*np.max(y_data)
                                        p = np.array([sigma, x_data[12:20][np.argmax(y_data[12:20])],-5,np.min(y_data),amp])
                                        try:
                                            popt, pcov = curve_fit(func, x_data, y_data, p)
                                            err = np.sum(np.abs(func(x_data, *popt) - y_data))
                                            results.append((err, popt))
                                            y_fit_data = func(x_data, *popt)
                                            chi_squared_pulsed = np.sum((long_median_total[23:-15] - y_fit_data)**2./err**2.)
                                            maximum_likelihood_pulsed = np.exp(-chi_squared_pulsed/2.)
                                            y_no_pulsed = np.ones(np.shape(y_data))*np.random.normal(np.mean(long_median_total), np.std(long_median_total))
                                            chi_squared_no_pulsed = np.sum((long_median_total[23:-15] - y_no_pulsed)**2./err**2.)
                                            maximum_likelihood_no_pulsed = np.exp(-chi_squared_no_pulsed/2.)
                                            pulsed_likelihoods.append(maximum_likelihood_pulsed)
                                            no_pulsed_likelihoods.append(maximum_likelihood_no_pulsed)
                                        except:
                                            pass
                                    bayes_factor = np.mean(pulsed_likelihoods)/np.mean(no_pulsed_likelihoods)
                                    print("bayes_factor="+str(bayes_factor)+" for $e$=["+str(e_bins[e_bin_it-1])+","+str(e_bins[e_bin_it])+"]")
                            
                                    err, popt = min(results, key=lambda x:x[0])
                                    y_fit= func(x, *popt)
                                except:
                                    popt=np.nan
                                #import pdb
                                #pdb.set_trace()
                            y_fits.append(popt)
                            #axes_dict[ax_label].plot(x, y_fit*(1.e4), 'k--')
                            '''
                            
                            if e_bin_it == len(e_bins) - 1:
                                #axes_dict[ax_label].legend(loc='center right', fontsize=args.text_font)#(loc='center left', bbox_to_anchor=(0.985, 0.5), fontsize=args.text_font)
                                #axes_dict[ax_label].legend(loc='center left', bbox_to_anchor=(0.985, 0.5), fontsize=args.text_font)
                                axes_dict[ax_label].legend(loc='center right', fontsize=args.text_font)
                            if plot_it == 2:
                                #fig.text(0.5, 0.07, "Orbital Phase ($\phi$)", va='center', ha='center')
                                plt.xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
                            plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
                            plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
                            #call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
                            e_bin_it = e_bin_it + 1
                            plot_it = plot_it + 1
                            
                            if args.pickle_dump != 'False':
                                #long_median_accretion = [np.array(long_median_accretion_1)*(1.e4), np.array(long_median_accretion_2)*(1.e4), np.array(long_median_total)*(1.e4)]
                                #yerr_tot = [np.array(yerr_1)*(1.e4), np.array(yerr_2)*(1.e4), np.array(yerr)*(1.e4)]
                                long_median_accretion = [np.array(long_median_accretion_1), np.array(long_median_accretion_2), np.array(long_median_total)]
                                yerr_tot = [np.array(yerr_1), np.array(yerr_2), np.array(yerr)]
                                component_accretion_pickle = file_name + '.pkl'
                                file = open(component_accretion_pickle, 'wb')
                                pickle.dump((phase_centers, long_median_accretion, yerr_tot, popt),file)
                                file.close()
                
                    n_lines = len(multiple_folds)
                    c_index = np.linspace(0.0, 0.95, n_lines)
                    alpha_values = np.linspace(1.0, 0.2, n_lines)
                    e_int = 1
                    plt.clf()
                    for i, shift in zip(c_index, multiple_folds):
                        try:
                            plt.plot(phase_centers, shift, color=plt.cm.gnuplot(i), label='$e$=['+str(e_bins[e_int-1])+','+str(e_bins[e_int])+']')#, ls='steps-mid')#, alpha=alpha_values[e_int])
                        except:
                            print("couldn't plot" + '$e$=['+str(e_bins[e_int-1])+','+str(e_bins[e_int])+']')
                        e_int = e_int + 1
                    plt.legend(loc='best')
                    plt.xlabel("Orbital Phase ($\phi$)")
                    #plt.ylabel("Accretion Rate ($10^{-4}$M$_\odot$/yr)")
                    plt.ylabel("Accretion Rate (M$_\odot$/yr)")
                    plt.xlim([0.0, 1.3])
                    #plt.ylim(bottom=0.0)
                    if args.pickle_dump != 'False':
                        folded_pickle = path + 'using_e_bins.pkl'
                        file = open(folded_pickle, 'wb')
                        pickle.dump((multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits),file)
                        file.close()
                    file_name = save_dir + 'multiple_folds_using_e_bins'
                    plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
                    plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
                break

