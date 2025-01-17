import numpy as np
import pickle
import yt
yt.enable_parallelism()
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
from mpi4py.MPI import COMM_WORLD as CW
import my_flash_fields as myf
import my_flash_module as mym
import pickle
import argparse
import os

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
#Ploting parameters
matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
max_time = [[10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000]]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

#---------------------------------------------------
#Define arguments
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-update", "--update_pickles", help="do you want to read the Flash output and update the pickles?", type=str, default='True')
    parser.add_argument("-lref", "--refinment_level", type=str, default='9')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
#input_dir = sys.argv[1]
args = parse_inputs()

#Mach_labels = ['0.0', '0.1', '0.2']
Mach_labels = ['0.0', '0.2']
Spin_labels = ['0.20', '0.25', '0.30', '0.35']

height_multiplier = 1.5
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Mach_labels), figsize=(single_col_width, height_multiplier*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymax = 0
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            for sink_id in sink_data.keys():
                #sink_id = list(sink_data.keys())[0]
                if np.isnan(form_time):
                    form_time = sink_data[sink_id]['time'][0]
                mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
                L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
                L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
                #L_tot = L_tot.in_units('kg*m**2/s')
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                    plot_time = time.in_units('yr')
                    plot_mass = mass.in_units('kg*m**2/s')
                else:
                    end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                    plot_time = time.in_units('yr')[:end_ind+1]
                    plot_L = L_tot.in_units('kg*m**2/s')[:end_ind+1]/1.e45
                if sink_id == list(sink_data.keys())[0]:
                    axs.flatten()[plot_it].plot(plot_time, plot_L, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
                else:
                    axs.flatten()[plot_it].plot(plot_time, plot_L, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)

            if len(sink_data.keys()) == 1:
                spin_up_percentage = 0
            else:
                prime_id = list(sink_data.keys())[0]
                second_id = list(sink_data.keys())[1]
                secondary_form_time = sink_data[second_id]['time'][0]
                secondary_form_ind = np.argmin(abs(sink_data[prime_id]['time'] - secondary_form_time))
                pre_sec_L = np.sqrt(sink_data[prime_id]['anglx'][secondary_form_ind-1]**2 + sink_data[prime_id]['angly'][secondary_form_ind-1]**2 + sink_data[prime_id]['anglz'][secondary_form_ind-1]**2)
                pre_sec_L_spec = pre_sec_L/sink_data[prime_id]['mass'][secondary_form_ind-1]
                post_sec_L = np.sqrt(sink_data[prime_id]['anglx'][-1]**2 + sink_data[prime_id]['angly'][-1]**2 + sink_data[prime_id]['anglz'][-1]**2)
                post_sec_L_spec = post_sec_L/sink_data[prime_id]['mass'][-1]
                DL = post_sec_L - pre_sec_L
                spin_up_percentage = DL/pre_sec_L * 100
                DL_spec = post_sec_L_spec - pre_sec_L_spec
                spin_up_percentage_spec = DL_spec/pre_sec_L_spec * 100
                print("Spin up: L=", spin_up_percentage, "%, L_spec=", spin_up_percentage_spec, "%")
            #Calculate pre spin adn post spin
        else:
            print("Couldn't open", single_pickle)
    
    axs.flatten()[plot_it].set_ylabel('L ($10^{45}kg\,m^2/s$)', labelpad=-0.2)
    
    if mach_lab == '0.0':
        mach_string = "No Turbulence \n($\mathcal{M}$="+mach_lab+")"
        mach_string_raw = r"{}".format(mach_string)
        time_text = axs.flatten()[plot_it].text(250, 0.85, mach_string_raw, va="center", ha="left", color='k', fontsize=font_size)
    else:
        mach_string = "With Turbulence \n($\mathcal{M}$="+mach_lab+")"
        mach_string_raw = r"{}".format(mach_string)
        time_text = axs.flatten()[plot_it].text(250, 1.2, mach_string_raw, va="center", ha="left", color='k', fontsize=font_size)

axs.flatten()[0].legend(loc='upper left', ncol=2)#, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].set_xlabel('Time ($yr$)', labelpad=-0.2)
yticklabels = axs.flatten()[1].get_yticklabels()
plt.setp(yticklabels[-1], visible=False)

axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim(bottom=0)

plt.savefig('Spin_init_spin_comp.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Mach_labels), figsize=(single_col_width, height_multiplier*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

spin_up_start = [[4500, 4500, 5250, 6750], [3750, 3500, 4250, 4250]]
spin_up_end = [[5400, 5450, 5750, 8000], [5750, 5300, 5500, 5750]]
peak_times = [[6288.53350698, 5337.46225949, 5961.44966664, 7499.56279945],
 [5684.37286105, 5479.48842751, 5771.10385454, 5758.7326983]]

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymax = 1.5222669979520504
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            plot_highlight = True
            for sink_id in sink_data.keys():
                #sink_id = list(sink_data.keys())[0]
                if np.isnan(form_time):
                    form_time = sink_data[sink_id]['time'][0]
                mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
                L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
                L_tot = yt.YTArray(L_tot/sink_data[sink_id]['mass'], 'cm**2/s')
                #L_tot = L_tot.in_units('m**2/s')
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                    plot_time = time.in_units('yr')
                    plot_mass = mass.in_units('m**2/s')
                else:
                    end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                    plot_time = time.in_units('yr')[:end_ind+1]
                    plot_L = L_tot.in_units('m**2/s')[:end_ind+1]/1.e15
                if sink_id == list(sink_data.keys())[0]:
                    axs.flatten()[plot_it].plot(plot_time, plot_L, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
                else:
                    axs.flatten()[plot_it].plot(plot_time, plot_L, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
                if plot_highlight == True:
                    highlight_start_time = spin_up_start[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
                    if np.isnan(highlight_start_time) == False:
                        highlight_end_time = spin_up_end[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
                        highlight_start_ind = np.argmin(abs(plot_time.value-highlight_start_time))
                        highlight_end_ind = np.argmin(abs(plot_time.value-highlight_end_time))
                        highlight_min = np.min(plot_L[highlight_start_ind:highlight_end_ind])
                        highlight_max = np.max(plot_L[highlight_start_ind:highlight_end_ind])
                        axs.flatten()[plot_it].axvspan(highlight_start_time, highlight_end_time, ymin=float(highlight_min.value)/ymax, ymax=float(highlight_max.value)/ymax, alpha=0.30, facecolor=colors[Spin_labels.index(spin_lab)])
                    plot_highlight = False
        else:
            print("Couldn't open", single_pickle)
    
    axs.flatten()[plot_it].set_ylabel('h ($10^{15}m^2/s$)', labelpad=-0.2)
    
    if mach_lab == '0.0':
        mach_string = "No Turbulence \n($\mathcal{M}$="+mach_lab+")"
        mach_string_raw = r"{}".format(mach_string)
        time_text = axs.flatten()[plot_it].text(250, 1.35, mach_string_raw, va="center", ha="left", color='k', fontsize=font_size)
    else:
        mach_string = "With Turbulence \n($\mathcal{M}$="+mach_lab+")"
        mach_string_raw = r"{}".format(mach_string)
        time_text = axs.flatten()[plot_it].text(250, 1.35, mach_string_raw, va="center", ha="left", color='k', fontsize=font_size)

axs.flatten()[0].legend(loc='lower right', ncol=2)#, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].set_xlabel('Time ($yr$)', labelpad=-0.2)
yticklabels = axs.flatten()[1].get_yticklabels()
plt.setp(yticklabels[-3], visible=False)

axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim(bottom = 0)

plt.savefig('Spin_init_spin_spec_comp.pdf', bbox_inches='tight', pad_inches=0.02)

#==========================================================================================================================
#Star stats
plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(two_col_width, 0.8*page_height), sharex=True, sharey='row')
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']

plot_quantity = ['mass', 'angular momentum', 'specific angular momentum', 'period']
radius = yt.YTQuantity(2, 'rsun')
M_eff = [0.7, 0.9]
L_eff = [0.01, 0.1]

for mach_lab in Mach_labels:
    axs.flatten()[6].axhline(y=2, color='k', linewidth=0.5)
    axs.flatten()[7].axhline(y=2, color='k', linewidth=0.5)
    for spin_lab in Spin_labels:
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9.pkl'
        file = open(single_pickle, 'rb')
        sink_data, line_counter = pickle.load(file)
        file.close()
        form_time = np.nan

        primary_ind = list(sink_data.keys())[0]
        form_time = sink_data[primary_ind]['time'][0]
        mass = yt.YTArray(sink_data[primary_ind]['mass'], 'g')
        L_tot = np.sqrt(sink_data[primary_ind]['anglx']**2 + sink_data[primary_ind]['angly']**2 + sink_data[primary_ind]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        time = sink_data[primary_ind]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - 10000))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        if mach_lab == '0.0':
            mach_string = "No Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[0].text(9500, 0.04, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
        else:
            mach_string = "With Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[1].text(9500, 0.04, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
            
        plot_it = Mach_labels.index(mach_lab) - 2
        for plot_q in plot_quantity:
            plot_it = plot_it + 2
            
            axs.flatten()[plot_it].tick_params(axis='x', direction='in', top=True)
            axs.flatten()[plot_it].tick_params(axis='y', direction='in', right=True)
            axs.flatten()[plot_it].minorticks_on()
            axs.flatten()[plot_it].tick_params(which='both', direction='in', axis='both', right=True, top=True)
            
            if plot_q == 'mass':
                m_star_lower = M_eff[0] * mass.in_units('g')[:end_ind+1]
                m_star_upper = M_eff[1] * mass.in_units('g')[:end_ind+1]
                m_star = ((m_star_lower + m_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, m_star.in_units('msun'), label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, m_star_lower.in_units('msun'), m_star_upper.in_units('msun'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$M_\star$ ($M_\odot$)', labelpad=-0.2)
                    axs.flatten()[plot_it].legend(loc='best', ncol=2, columnspacing=0.8)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(m_star_upper.in_units('msun'))])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'angular momentum':
                l_star_lower = L_eff[0] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star_upper = L_eff[1] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star = ((l_star_lower + l_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, l_star.in_units('kg*m**2/s')/1e44, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, l_star_lower.in_units('kg*m**2/s')/1e44, l_star_upper.in_units('kg*m**2/s')/1e44, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$L_\star$ ($10^{44}kg\,m^2/s$)', labelpad=-0.2)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(l_star_upper.in_units('kg*m**2/s')/1e44)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'specific angular momentum':
                h_star_lower = l_star_lower/m_star_upper
                h_star_upper = l_star_upper/m_star_lower
                h_star = ((h_star_lower + h_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, h_star.in_units('m**2/s')/1e14, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, h_star_lower.in_units('m**2/s')/1e14, h_star_upper.in_units('m**2/s')/1e14, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$h_\star$ ($10^{14}m^2/s$)', labelpad=-0.2)
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-2], visible=False)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(h_star_upper.in_units('m**2/s')/1e14)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'period':
                P_star_lower = ((4*np.pi)/5) * (radius.in_units('m')**2)/h_star_upper
                P_star_upper = ((4*np.pi)/5) * (radius.in_units('m')**2)/h_star_lower
                P_star = ((P_star_lower + P_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, P_star.in_units('day'), linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, P_star_lower.in_units('day'), P_star_upper.in_units('day'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].set_xlabel('Time (yr)', labelpad=-0.2)
                axs.flatten()[plot_it].set_ylim([0, 5])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$P_\star$ (days)')
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-1], visible=False)

axs.flatten()[plot_it-1].set_xlim([0, 10000])

plt.savefig('period.pdf', bbox_inches='tight', pad_inches=0.02)

#=======================================================================================
#From resolution study

#Star stats
plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(two_col_width, 0.8*page_height), sharex=True, sharey='row')
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']

plot_quantity = ['mass', 'angular momentum', 'specific angular momentum', 'period']
radius = yt.YTQuantity(2, 'rsun')

M_eff = [0.05210517772546909, 0.20648041032355596]
L_eff = [0.0014253666342484849, 0.005551756458477013]
h_eff = [0.014254361678943728, 0.025101068720250228]
T_eff = [0.0001263524949463752, 0.00020639462094946843]

for mach_lab in Mach_labels:
    axs.flatten()[6].axhline(y=2, color='k', linewidth=0.5)
    axs.flatten()[7].axhline(y=2, color='k', linewidth=0.5)
    for spin_lab in Spin_labels:
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9.pkl'
        file = open(single_pickle, 'rb')
        sink_data, line_counter = pickle.load(file)
        file.close()
        form_time = np.nan

        primary_ind = list(sink_data.keys())[0]
        form_time = sink_data[primary_ind]['time'][0]
        mass = yt.YTArray(sink_data[primary_ind]['mass'], 'g')
        L_tot = np.sqrt(sink_data[primary_ind]['anglx']**2 + sink_data[primary_ind]['angly']**2 + sink_data[primary_ind]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        time = sink_data[primary_ind]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - 10000))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        if mach_lab == '0.0':
            mach_string = "No Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[0].text(9500, 0.04, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
        else:
            mach_string = "With Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[1].text(9500, 0.04, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
            
        plot_it = Mach_labels.index(mach_lab) - 2
        for plot_q in plot_quantity:
            plot_it = plot_it + 2
            
            axs.flatten()[plot_it].tick_params(axis='x', direction='in', top=True)
            axs.flatten()[plot_it].tick_params(axis='y', direction='in', right=True)
            axs.flatten()[plot_it].minorticks_on()
            axs.flatten()[plot_it].tick_params(which='both', direction='in', axis='both', right=True, top=True)
            
            if plot_q == 'mass':
                m_star_lower = M_eff[0] * mass.in_units('g')[:end_ind+1]
                m_star_upper = M_eff[1] * mass.in_units('g')[:end_ind+1]
                m_star = ((m_star_lower + m_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, m_star.in_units('msun'), label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, m_star_lower.in_units('msun'), m_star_upper.in_units('msun'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$M_\star$ ($M_\odot$)', labelpad=-0.2)
                    axs.flatten()[plot_it].legend(loc='best', ncol=2, columnspacing=0.8)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(m_star_upper.in_units('msun'))])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'angular momentum':
                l_star_lower = L_eff[0] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star_upper = L_eff[1] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star = ((l_star_lower + l_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, l_star.in_units('kg*m**2/s')/1e44, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, l_star_lower.in_units('kg*m**2/s')/1e44, l_star_upper.in_units('kg*m**2/s')/1e44, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$L_\star$ ($10^{44}kg\,m^2/s$)', labelpad=-0.2)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(l_star_upper.in_units('kg*m**2/s')/1e44)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'specific angular momentum':
                h_star_lower = l_star_lower/m_star_upper
                h_star_upper = l_star_upper/m_star_lower
                h_star = ((h_star_lower + h_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, h_star.in_units('m**2/s')/1e14, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, h_star_lower.in_units('m**2/s')/1e14, h_star_upper.in_units('m**2/s')/1e14, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$h_\star$ ($10^{14}m^2/s$)', labelpad=-0.2)
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-2], visible=False)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(h_star_upper.in_units('m**2/s')/1e14)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'period':
                P_star_lower = ((4*np.pi)/5) * (radius.in_units('m')**2)/h_star_upper
                P_star_upper = ((4*np.pi)/5) * (radius.in_units('m')**2)/h_star_lower
                P_star = ((P_star_lower + P_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, P_star.in_units('day'), linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, P_star_lower.in_units('day'), P_star_upper.in_units('day'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].set_xlabel('Time (yr)', labelpad=-0.2)
                axs.flatten()[plot_it].set_ylim([0, 5])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$P_\star$ (days)')
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-1], visible=False)

axs.flatten()[plot_it-1].set_xlim([0, 10000])

plt.savefig('period_resolution_study_extrapolation.pdf', bbox_inches='tight', pad_inches=0.02)

#=======================================================================================
#From resolution study

#Star stats
plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(two_col_width, 0.8*page_height), sharex=True, sharey='row')
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']

plot_quantity = ['mass', 'angular momentum', 'specific angular momentum', 'period']
radius = yt.YTQuantity(2, 'rsun')
r_sink_9 = yt.YTQuantity(4.89593797, 'AU')

M_eff = [0.05210517772546909, 0.20648041032355596]
L_eff = [0.0014253666342484849, 0.005551756458477013]
h_eff = [0.014254361678943728, 0.025101068720250228]
T_eff = [0.0001263524949463752, 0.00020639462094946843]

for mach_lab in Mach_labels:
    axs.flatten()[6].axhline(y=2, color='k', linewidth=0.5)
    axs.flatten()[7].axhline(y=2, color='k', linewidth=0.5)
    for spin_lab in Spin_labels:
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9.pkl'
        file = open(single_pickle, 'rb')
        sink_data, line_counter = pickle.load(file)
        file.close()
        form_time = np.nan

        primary_ind = list(sink_data.keys())[0]
        form_time = sink_data[primary_ind]['time'][0]
        mass = yt.YTArray(sink_data[primary_ind]['mass'], 'g')
        L_tot = np.sqrt(sink_data[primary_ind]['anglx']**2 + sink_data[primary_ind]['angly']**2 + sink_data[primary_ind]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        h_tot = L_tot.in_units('g*m**2/s')/mass.in_units('g')
        T_sink = ((4*np.pi)/5) * (r_sink_9.in_units('m')**2)/h_tot
        time = sink_data[primary_ind]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - 10000))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        if mach_lab == '0.0':
            mach_string = "No Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[0].text(9750, 0.005, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
        else:
            mach_string = "With Turbulence ($\mathcal{M}$="+mach_lab+")"
            mach_string_raw = r"{}".format(mach_string)
            time_text = axs.flatten()[1].text(9750, 0.005, mach_string_raw, va="center", ha="right", color='k', fontsize=font_size)
            
        plot_it = Mach_labels.index(mach_lab) - 2
        for plot_q in plot_quantity:
            plot_it = plot_it + 2
            
            axs.flatten()[plot_it].tick_params(axis='x', direction='in', top=True)
            axs.flatten()[plot_it].tick_params(axis='y', direction='in', right=True)
            axs.flatten()[plot_it].minorticks_on()
            axs.flatten()[plot_it].tick_params(which='both', direction='in', axis='both', right=True, top=True)
            
            if plot_q == 'mass':
                m_star_lower = M_eff[0] * mass.in_units('g')[:end_ind+1]
                m_star_upper = M_eff[1] * mass.in_units('g')[:end_ind+1]
                m_star = ((m_star_lower + m_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, m_star.in_units('msun'), label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, m_star_lower.in_units('msun'), m_star_upper.in_units('msun'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$M_\star$ ($M_\odot$)', labelpad=-0.2)
                    axs.flatten()[plot_it].legend(loc='upper left', ncol=2, columnspacing=0.8)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(m_star_upper.in_units('msun'))])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'angular momentum':
                l_star_lower = L_eff[0] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star_upper = L_eff[1] * L_tot.in_units('g*cm**2/s')[:end_ind+1]
                l_star = ((l_star_lower + l_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, l_star.in_units('kg*m**2/s')/1e44, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, l_star_lower.in_units('kg*m**2/s')/1e44, l_star_upper.in_units('kg*m**2/s')/1e44, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$L_\star$ ($10^{44}kg\,m^2/s$)', labelpad=-0.2)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(l_star_upper.in_units('kg*m**2/s')/1e44)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'specific angular momentum':
                h_star_lower = h_eff[0] * h_tot.in_units('m**2/s')[:end_ind+1]
                h_star_upper = h_eff[1] * h_tot.in_units('m**2/s')[:end_ind+1]
                h_star = ((h_star_lower + h_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, h_star.in_units('m**2/s')/1e14, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, h_star_lower.in_units('m**2/s')/1e14, h_star_upper.in_units('m**2/s')/1e14, alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$h_\star$ ($10^{14}m^2/s$)', labelpad=-0.2)
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-2], visible=False)
                if mach_lab == '0.2' and spin_lab == '0.20':
                    axs.flatten()[plot_it].set_ylim([0, np.max(h_star_upper.in_units('m**2/s')/1e14)])
                    axs.flatten()[plot_it].set_ylim(bottom=0)
            if plot_q == 'period':
                P_star_lower = T_eff[0] * T_sink.in_units('day')[:end_ind+1]
                P_star_upper = T_eff[1] * T_sink.in_units('day')[:end_ind+1]
                P_star = ((P_star_lower + P_star_upper)/2)
                axs.flatten()[plot_it].plot(plot_time, P_star.in_units('day'), linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].fill_between(plot_time, P_star_lower.in_units('day'), P_star_upper.in_units('day'), alpha=0.2, color=colors[Spin_labels.index(spin_lab)])
                axs.flatten()[plot_it].set_xlabel('Time (yr)', labelpad=-0.2)
                axs.flatten()[plot_it].set_ylim([0, 5])
                if mach_lab == '0.0':
                    axs.flatten()[plot_it].set_ylabel('$P_\star$ (days)')
                    if spin_lab == '0.35':
                        yticklabels = axs.flatten()[plot_it].get_yticklabels()
                        plt.setp(yticklabels[-1], visible=False)

axs.flatten()[plot_it-1].set_xlim([0, 10000])

plt.savefig('period_resolution_study.pdf', bbox_inches='tight', pad_inches=0.02)
            
