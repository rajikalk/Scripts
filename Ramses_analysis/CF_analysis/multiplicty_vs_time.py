import numpy as np
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW
import collections

f_acc= 0.5

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-sim_G", "--simulation_G", type=str, default='')
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("-start_ind", "--starting_ind", help="Do you want to start the analysis at a particular starting ind?", type=int, default=None)
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=55.29)
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.07)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

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

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

datadir = sys.argv[1]
savedir = sys.argv[2]
#=====================================================================================================
#Create units override

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_G == '':
    simulation_density_id = args.global_data_pickle_file.split('/G')[-1].split('/')[0]
else:
    simulation_density_id = args.simulation_G

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
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

#====================================================================================================================================================
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

sys.stdout.flush()
CW.Barrier()

if args.make_plots_only == 'False':
    file_open = open(args.global_data_pickle_file, 'rb')
    try:
        global_data = pickle.load(file_open,encoding="latin1")
    except:
        file_open.close()
        import pickle5 as pickle
        file_open = open(args.global_data_pickle_file, 'rb')
        global_data = pickle.load(file_open,encoding="latin1")
    file_open.close()
    print('Loaded global pickle data')
        
    SFE = np.sum(global_data['m'], axis=1)
    SFE_vals = [0.01, 0.02, 0.03, 0.04, 0.05]
    MF_plot = []
    MF_plot_err = []
    SFE_window = 0.001
    time_its = []
    mass_bins = np.logspace(-1.5, 1.5, 10)
    for SFE_val in SFE_vals:
        start_SFE = SFE_val - SFE_window
        end_SFE = SFE_val + SFE_window
        start_time_it = np.argmin(abs(SFE-start_SFE))
        end_time_it = np.argmin(abs(SFE-end_SFE))
        time_its = np.arange(start_time_it, end_time_it)
        
        sys.stdout.flush()
        CW.Barrier()
        
        MF_arrays = []
        rit = -1
        for time_it in time_its:
            rit = rit + 1
            if rit == size:
                rit = 0
            if rank == rit:
                n_stars = np.where(global_data['m'][time_it]>0)[0]
                abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                mass = np.array(global_data['m'][time_it][n_stars])
                time = global_data['time'][time_it][n_stars][0]
                
                #True multiplicity
                S = pr.Sink()
                S._jet_factor = 1.
                S._scale_l = scale_l.value
                S._scale_v = scale_v.value
                S._scale_t = scale_t.value
                S._scale_d = scale_d.value
                S._time = yt.YTArray(time, '')
                S._abspos = yt.YTArray(abspos, '')
                S._absvel = yt.YTArray(absvel, '')
                S._mass = yt.YTArray(mass, '')
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, Grho=Grho, max_iter=100)
                top_inds = np.where(res['topSystem'])[0]
                
                IMF = np.histogram(units_override['mass_unit'][0]*mass, bins=np.logspace(-1.5, 1.5, 10))[0]
                primary_masses = []
                for top_ind in top_inds:
                    top_sys_ids = flatten(losi(top_ind, res))
                    primary_masses.append(np.max(res['mass'][top_sys_ids]))
                primary_hist = np.histogram(primary_masses, bins=np.logspace(-1.5, 1.5, 10))[0]
                MF_arrays.append((primary_hist/IMF))
                print('found MF for time_it', time_it-time_its[0], 'of', time_its[-1]-time_its[0], 'on rank', rank)
                
        sys.stdout.flush()
        CW.Barrier()
        
        file = open('MF_rank_'+str(rank)+'.pkl', 'wb')
        pickle.dump((MF_arrays), file)
        file.close()
        
        sys.stdout.flush()
        CW.Barrier()
        
        if rank == 0:
            MF_pickles = glob.glob('MF_rank_*.pkl')
            MF_all = []
            for MF_pickle in MF_pickles:
                file = open(MF_pickle, 'rb')
                MF = pickle.load(file)
                file.close()
                MF_all = MF_all + MF
                os.remove(MF_pickle)
            
            MF_median = np.nanmedian(MF_all, axis=0)
            MF_mean = np.nanmean(MF_all, axis=0)
            MF_std = np.nanstd(MF_all, axis=0)
            MF_err = [MF_median-(MF_mean-MF_std), (MF_mean+MF_std)-MF_median]
            MF_plot.append(MF_median)
            MF_plot_err.append(MF_err)
        
        sys.stdout.flush()
        CW.Barrier()

    if rank == 0:
        file = open('MF_plot.pkl', 'wb')
        pickle.dump((SFE_vals, MF_plot, MF_plot_err), file)
        file.close()
    
sys.stdout.flush()
CW.Barrier()

if rank == 0:
    file = open('MF_plot.pkl', 'rb')
    SFE_vals, MF_plot, MF_plot_err = pickle.load(file)
    file.close()

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
    #matplotlib.rcParams['text.latex.preamble'] = [
    #       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    #       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    #       r'\usepackage{helvet}',    # set the normal font here
    #       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    #       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    #]
    
    import matplotlib.pylab as pl
    colors = pl.cm.cool(np.linspace(0,1,len(SFE_vals)))
    
    import matplotlib.pyplot as plt

    mass_bins = np.logspace(-1.5, 1.5, 10)
    bin_centres = np.log10((mass_bins[1:] + mass_bins[:-1])/2)

    plt.clf()
    plt.figure(figsize=(single_col_width, single_col_width))
    for SFE_it in range(len(SFE_vals)):
        plt.errorbar(bin_centres, MF_plot[SFE_it], np.array(MF_plot_err[SFE_it]), color=colors[SFE_it], label='SFE='+str(SFE_vals[SFE_it]))
    plt.legend(loc='best')
    plt.xlabel('Log$_{10}$ Mass')
    plt.ylabel('Multiplicity Fraction')
    #plt.ylim(bottom=0.0)
    plt.tick_params(which='both', direction='in')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=10)
    plt.savefig('MF_evolution.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)
    print('made MF_evolution.pdf')
