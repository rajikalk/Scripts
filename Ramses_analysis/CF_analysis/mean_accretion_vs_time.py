import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
from mpi4py.MPI import COMM_WORLD as CW
import sys
import collections
import os
import pickle
import gc
import matplotlib

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

f_acc= 0.5

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def accretion(sink_inds, time_ind):
    """
    Calculates the accretion of the given indeexes
    """
    global Accretion_array
    M_dot = Accretion_array[time_ind, sink_inds]
    return M_dot
    
def luminosity(global_data, sink_inds, global_ind):
    """
    Calculates the luminosity of the given indexes
    """
    global f_acc
    radius = yt.YTQuantity(2.0, 'rsun')
    M_dot = accretion(sink_inds, global_ind)
    M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    return L_tot

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-update", "--update_pickles", help="Do you want to remake the pickles?", type=str, default='True')
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=55.29)
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.09)
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-max_iter", "--max_iterations", help="How many iterations for multiplicity analysis", type=int, default=20)
    parser.add_argument("-start_time_it", "--start_time_index", help="What time index do you want to start at? mostly for diagnostic reasons.", type=int, default=0)
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifeimteimte threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-plot_super", "--plot_superplot", help="do you want to plot the superplot?", type=str, default='True')
    parser.add_argument("-use_com_sep", "--use_separation_of_the_coms", help="Do you want to same the separation as the separation of the center of masses instead of relative to a particular sink?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

global_data_pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G50/stars_imf_G50.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G100/256/stars_imf_G100.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G125/stars_imf_G125.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G150/stars_imf_G150.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G200/stars_imf_G200.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G400/stars_imf_G400.pkl"]

labels = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
line_styles = [':', (0, (3, 1, 1, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), '--', '-']

plt.clf()
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.7*single_col_width))
pit = -1
for global_data_pickle_file in global_data_pickle_files:
    print("reading file", global_data_pickle_file)
    pit = pit+1
    #Set units
    units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

    simulation_density_id = global_data_pickle_file.split('/')[-1].split('_')[0][1:]
    try:
        simulation_density_id_int = int(simulation_density_id)
    except:
        simulation_density_id = global_data_pickle_file.split('_G')[-1].split('.pkl')[0]

    if simulation_density_id == '50':
        Grho=50
        units_override.update({"mass_unit":(1500,"Msun")})
    elif simulation_density_id == '100':
        Grho=100
        units_override.update({"mass_unit":(3000,"Msun")})
    elif simulation_density_id == '125':
        Grho=125
        units_override.update({"mass_unit":(3750,"Msun")})
    elif simulation_density_id == '150':
        Grho=150
        units_override.update({"mass_unit":(4500,"Msun")})
    elif simulation_density_id == '200':
        Grho=200
        units_override.update({"mass_unit":(6000,"Msun")})
    elif simulation_density_id == '400':
        Grho=400
        units_override.update({"mass_unit":(12000,"Msun")})
    else:
        print("MASS UNIT NOT SET")
        import pdb
        pdb.set_trace()
        

    units_override.update({"density_unit":(units_override['mass_unit'][0]/(units_override['length_unit'][0]**3), "Msun/pc**3")})
        
    scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
    scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
    scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
    scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

    units={}
    for key in units_override.keys():
        units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
    del units_override
    gc.collect()

    sys.stdout.flush()
    CW.Barrier()

    #loading global data
    file_open = open(global_data_pickle_file, 'rb')
    try:
        global_data = pickle.load(file_open,encoding="latin1")
    except:
        file_open.close()
        import pickle5 as pickle
        file_open = open(global_data_pickle_file, 'rb')
        global_data = pickle.load(file_open,encoding="latin1")
    file_open.close()
    dm = global_data['dm']*units['mass_unit'].in_units('Msun')
    dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')
    Accretion_array = dm/dt
    print('loaded global data')

    #Plot mean TOTAL accretion rates over time
    Mean_acc = np.nanmean(Accretion_array, axis=1)
    SFE_arr = np.nansum(global_data['m'], axis=1)
    plt.semilogy(SFE_arr, Mean_acc, label=labels[pit], color=colors[pit], linestyle=line_styles[pit])
plt.legend(ncol=2, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
plt.xlim([0, 0.05])
plt.xlabel('SFE')
plt.ylabel("Accretion rate (Msun/yr)")
savename = "mean_acc.pdf"
plt.savefig(savename)
