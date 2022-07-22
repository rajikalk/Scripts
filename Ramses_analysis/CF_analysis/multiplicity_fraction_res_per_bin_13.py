import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import collections
import matplotlib as mpl
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW
import matplotlib.gridspec as gridspec
from scipy.stats import norm
from scipy.optimize import curve_fit

#Define globals
f_acc= 0.5
Accretion_array = []

Perseus_L_lims = [0.07, 55.29]
Orion_L_lims = [0.02, 1404.97]

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-proj", "--projected_separation", help="do you want to use projected separation instead of true separation?", default="False", type=str)
    parser.add_argument("-ax", "--axis", help="what axis do you want to project separation onto?", default='z', type=str)
    parser.add_argument("-preffix", "--figure_prefix", help="Do you want to give saved figures a preffix?", default="", type=str)
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-verbose", "--verbose_printing", help="Would you like to print debug lines?", type=str, default='False')
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=55.29)#120?
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.07)
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-use_midpoint", "--use_midpoint_separation", help="Do you want to use the midpoint separation instread of separation", type=str, default='')
    parser.add_argument("-lifetime", "--lifetime_threshold", help="What life time threshold do you want to consider when making Luminosity histogram", type=float, default=10000)
    parser.add_argument("-proj_vec", "--projection_vector", help="What projection vector do you want to use?", type=str, default='')
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=10000)
    parser.add_argument("-use_t_s", "--use_t_spread", help="Do you want to define the time spread using t_spread instead of a spread defined by the method?", type=str, default='True')
    parser.add_argument("-match_meth", "--match_method", help="How do you want to select times? 1. SFE, 2.SPE/t_ff, 3. No. visible stars, 4. Total Accreted Mass", type=int, default=1)
    parser.add_argument("-start_ind", "--starting_ind", help="Do you want to start the analysis at a particular starting ind?", type=int, default=None)
    parser.add_argument("-n_vis_thres", "--n_visible_threshold", help="what threshold do you want to use for number fo stars?", type=int, default=106)
    parser.add_argument("-SFE_thres", "--SFE_threshold", help="What threshold do you want to use for the SFE?", type=float, default=0.045)
    parser.add_argument("-SFE_n_thres", "--SFE_n_threshold", help="What threshod do you want to use for the SFE_n (defined by ARce et al, 2010)", type=float, default=0.022)
    parser.add_argument("-thres_spread", "--threshold_spread", help="Over what spread (as a fraction of the threshold) do you want to intergrate?", type=float, default=0.1)
    parser.add_argument("-entire_sim", "--integrate_over_entire_sim", help="Do you want to integrate over the entire sim?", type=str, default="False")
    parser.add_argument("-sim_G", "--simulation_G", type=str, default='')
    parser.add_argument("-debug", "--debugging", help="This flag is to stop at PDB steps", type=str, default="False")
    parser.add_argument("-vis_only", "--visible_only", help="Do you only want to feed the visible stars into multiplicity analysis?", type=str, default="False")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
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
        
def is_hierarchical(sys_structure):
    close_braket = False
    for char in str(sys_structure):
        if char == ']':
            close_braket = True
        if char == '[' and close_braket == True:
            is_hierarchical = False
            break
        else:
            is_hierarchical = True
    return is_hierarchical

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
    

#=====================================================================================================

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

datadir = sys.argv[1]
savedir = sys.argv[2]

#=====================================================================================================
#Tobin data
L_bins = np.logspace(-1.25,3.5,20)
S_bins = np.logspace(0.75,4,14)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
#=====================================================================================================
#Create units override

if len(args.projection_vector) > 0:
    proj_vector = np.array(eval(args.projection_vector))
    proj_length = np.sqrt(np.sum(proj_vector**2))
    proj_unit = proj_vector/proj_length
else:
    proj_unit = []

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

sys.stdout.flush()
CW.Barrier()

#What do I need to save?
Times = []
SFE = []
N_sys_total = []
CF_arrays = []
Sink_Luminosities = {}
Sink_Accretion = {}

#Overwrite saved variables if need be
pickle_file = savedir + args.pickled_file
if rank == 0:
    pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
    if len(pickle_files) > 0:
        Times_full = []
        SFE_full = []
        CF_per_bin_full = []
        n_systems_full = []
        Sink_Luminosities_full = {}
        Sink_Accretion_full = {}
        for pick_file in pickle_files:
            file = open(pick_file, 'rb')
            CF_arrays, N_sys_total, Times, SFE, Sink_Luminosities, Sink_Accretion = pickle.load(file)
            file.close()
            for L_key in Sink_Luminosities.keys():
                if L_key not in Sink_Luminosities_full.keys():
                    Sink_Luminosities_full.update({L_key:Sink_Luminosities[L_key]})
                    Sink_Accretion_full.update({L_key:Sink_Accretion[L_key]})
                else:
                    Sink_Luminosities_full[L_key] = Sink_Luminosities_full[L_key] + Sink_Luminosities[L_key]
                    Sink_Accretion_full[L_key] = Sink_Accretion_full[L_key] + Sink_Accretion[L_key]
            Times_full = Times_full + Times
            SFE_full = SFE_full + SFE
            CF_per_bin_full = CF_per_bin_full + CF_arrays
            n_systems_full = n_systems_full + N_sys_total
        
        #Let's sort the data
        for L_key in Sink_Luminosities_full:
            sorted_array = np.array(Sink_Luminosities_full[L_key])
            dict_sort = np.argsort(sorted_array.T[0])
            sorted_array.T[0] = np.array(Sink_Luminosities_full[L_key]).T[0][dict_sort]
            sorted_array.T[1] = np.array(Sink_Luminosities_full[L_key]).T[1][dict_sort]
            Sink_Luminosities[L_key] = sorted_array.tolist()
            sorted_array.T[0] = np.array(Sink_Accretion_full[L_key]).T[0][dict_sort]
            sorted_array.T[1] = np.array(Sink_Accretion_full[L_key]).T[1][dict_sort]
            Sink_Accretion[L_key] = sorted_array.tolist()
        sorted_inds = np.argsort(Times_full)
        Times = np.array(Times_full)[sorted_inds].tolist()
        SFE = np.array(SFE_full)[sorted_inds].tolist()
        CF_arrays = np.array(CF_per_bin_full)[sorted_inds].tolist()
        N_sys_total = np.array(n_systems_full)[sorted_inds].tolist()
    else:
        try:
            file = open(pickle_file+'.pkl', 'rb')
            Times, SFE, CF_arrays, N_sys_total, Sink_Luminosities_full, Sink_Accretion_full = pickle.load(file)
            file.close()
        except:
            print("No pickles seem to exist, so starting out fresh")

sys.stdout.flush()
CW.Barrier()

#==========================================================================================
#Analysis parameters: Do you want to check for boundedness, or have luminosity limits?
if args.bound_check == 'True':
    bound_check = True
else:
    bound_check = False
    
if args.projected_separation == 'True':
    multiplicity_analysis_projection = True
    use_mid_point_sep = True
else:
    multiplicity_analysis_projection = False
    use_mid_point_sep = False

if args.use_midpoint_separation != '':
    if args.use_midpoint_separation == 'True':
        use_mid_point_sep = True
    else:
        use_mid_point_sep = False

if use_mid_point_sep:
    sep_key = 'midpointSep'
else:
    sep_key = 'separation'

#Calculate variables
luminosity_lower_limit = args.lower_L_limit# 0.04 #0.01
#luminosity_upper_limit = args.upper_L_limit
accretion_limit = args.accretion_limit

#==============================================================================================================================

file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
'''
dm = (global_data['m'][2:] - global_data['m'][:-2])*units['mass_unit'].in_units('Msun')
dt = (global_data['time'][2:] - global_data['time'][:-2])*units['time_unit'].in_units('yr')
Accretion_array = np.hstack((np.nan*np.zeros(((dm/dt).shape[1], 1)), (dm/dt).T, np.nan*np.zeros(((dm/dt).shape[1], 1))))
Accretion_array = yt.YTArray(Accretion_array.T, 'Msun/yr')
'''
dm = global_data['dm']*units['mass_unit'].in_units('Msun')
dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')
Accretion_array = dm/dt

print('Loaded global pickle data')

#dt == integration window, if you don't want to integrate over the entire simulation
dt = yt.YTQuantity(args.time_spread, 'yr')

if args.integrate_over_entire_sim == "True":
    time_bounds = [global_data['time'].T[0][0]*units['time_unit'].in_units('yr'), global_data['time'].T[0][-1]*units['time_unit'].in_units('yr')]
elif args.match_method == 1:
    #Match at SFE of 4.9%
    SFE_value = args.SFE_threshold
    SFE = (np.sum(global_data['m'], axis=1)*units['mass_unit'].value)/units['mass_unit'].value
    if args.use_t_spread == 'True':
        SFE_ind = np.argmin(np.abs(SFE_value-SFE))
        SFE_t_ff = global_data['time'].T[0][SFE_ind]*units['time_unit'].in_units('yr')
        time_bounds = [(SFE_t_ff-dt),(SFE_t_ff+dt)]
    else:
        SFE_spread = SFE_value*args.threshold_spread
        SFE_min_ind = np.argmin(np.abs((SFE_value-SFE_spread)-SFE))
        SFE_t_ff_min = global_data['time'].T[0][SFE_min_ind]*units['time_unit'].in_units('yr')
        SFE_max_ind = np.argmin(np.abs((SFE_value+SFE_spread)-SFE))
        SFE_t_ff_max = global_data['time'].T[0][SFE_max_ind]*units['time_unit'].in_units('yr')
        time_bounds = [SFE_t_ff_min*units['time_unit'].in_units('yr'),SFE_t_ff_max*units['time_unit'].in_units('yr')]
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), SFE)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("SFE")
        plt.axhline(y=SFE_value)
        plt.axvline(x=time_bounds[0])
        plt.axvline(x=time_bounds[1])
        plt.savefig("SFE.jpg")
elif args.match_method == 2:
    #Match at SFE of 3.4%
    SFE_n_value = args.SFE_n_threshold
    SFE = (np.sum(global_data['m'], axis=1)*units['mass_unit'].value)/units['mass_unit'].value
    SFE_n = SFE/(global_data['time'].T[0]*units['time_unit'].in_units('yr').value)
    if args.use_t_spread == 'True':
        SFE_n_ind = np.argmin(np.abs(SFE_n_value - SFE_n))
        SFE_t_ff = global_data['time'].T[0][SFE_n_ind]*units['time_unit'].in_units('yr')
        time_bounds = [(SFE_t_ff-dt),(SFE_t_ff+dt)]
    else:
        SFE_n_spread = SFE_n_value*args.threshold_spread
        SFE_n_min = np.argmin(np.abs((SFE_n_value-SFE_n_spread) - SFE_n))
        SFE_n_t_ff_min = global_data['time'].T[0][SFE_n_min]*units['time_unit'].in_units('yr')
        SFE_n_max = np.argmin(np.abs((SFE_n_value+SFE_n_spread) - SFE_n))
        SFE_n_t_ff_max = global_data['time'].T[0][SFE_n_max]*units['time_unit'].in_units('yr')
        time_bounds = [SFE_n_t_ff_min*units['time_unit'].in_units('yr'),SFE_n_t_ff_max*units['time_unit'].in_units('yr')]
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), SFE_n)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("SFE$_n$")
        plt.axhline(y=SFE_n_value)
        plt.axvline(x=time_bounds[0])
        plt.axvline(x=time_bounds[1])
        plt.savefig("SFE_n.jpg")
elif args.match_method == 3:
    N_vis_val = args.n_visible_threshold#115#70#115
    N_prev = 0
    time_yrs = (global_data['time'].T[0] - global_data['time'].T[0][0])*units['time_unit'].in_units('yr')
    total_stars = np.sum(global_data['m']>0, axis=1)
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), total_stars)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("N$_{stars}$")
        plt.savefig("total_stars.jpg")
    if args.use_t_spread == 'True':
        potential_inds = []
        N_vis_array = []
        for time in global_data['time'].T[0]:
            time_yt = yt.YTArray(time*scale_t, 's')
            global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
            sink_inds = np.where(global_data['m'][global_ind]>0)[0]
            L_tot = luminosity(global_data, sink_inds, global_ind)
            M_dot = accretion(sink_inds, global_ind)
            vis_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
            N_vis = len(vis_inds)
            N_vis_array.append(N_vis)
            if N_vis == N_vis_val:
                potential_inds.append(global_ind)
            elif N_vis_val > N_prev and N_vis_val < N_vis:
                potential_inds.append(global_ind)
            N_prev = N_vis
        if rank == 0:
            plt.clf()
            plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), total_stars)
            plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), N_vis_array)
            plt.xlabel("Time (t$_{ff}$)")
            plt.ylabel("N$_{stars}$")
            plt.axhline(y=N_vis_val)
            plt.savefig("visible_stars.jpg")
        N_ind = int(len(potential_inds)/2)
        
        N_t_ff = global_data['time'].T[0][potential_inds[N_ind]]
        time_bounds = [(N_t_ff-dt)*units['time_unit'].in_units('yr'),(N_t_ff+dt)*units['time_unit'].in_units('yr')]
    else:
        N_vis_spread = N_vis_val*args.threshold_spread
        potential_min_inds = []
        potential_max_inds = []
        N_vis_array = []
        for time in global_data['time'].T[0]:
            time_yt = yt.YTArray(time*scale_t, 's')
            global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
            sink_inds = np.where(global_data['m'][global_ind]>0)[0]
            L_tot = luminosity(global_data, sink_inds, global_ind)
            M_dot = accretion(sink_inds, global_ind)
            vis_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
            N_vis = len(vis_inds)
            N_vis_array.append(N_vis)
            if N_vis == (N_vis_val-N_vis_spread):
                potential_min_inds.append(global_ind)
            elif (N_vis_val-N_vis_spread)> N_prev and (N_vis_val-N_vis_spread)<N_vis:
                potential_min_inds.append(global_ind)
            if N_vis == (N_vis_val+N_vis_spread):
                potential_max_inds.append(global_ind)
            elif (N_vis_val+N_vis_spread) > N_prev and (N_vis_val+N_vis_spread)<N_vis:
                potential_max_inds.append(global_ind)
            N_prev = N_vis
        if rank == 0:
            plt.clf()
            plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), total_stars)
            plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), N_vis_array)
            plt.xlabel("Time (t$_{ff}$)")
            plt.ylabel("N$_{stars}$")
            plt.axhline(y=N_vis_val)
            plt.savefig("visible_stars.jpg")
        N_min = int(len(potential_min_inds)/2)
        N_max = int(len(potential_max_inds)/2)
        N_t_ff_min = global_data['time'].T[0][potential_min_inds[N_min]]
        try:
            N_t_ff_max = global_data['time'].T[0][potential_max_inds[N_max]]
        except:
            N_t_ff_max = global_data['time'].T[0][-1]
        time_bounds = [N_t_ff_min*units['time_unit'].in_units('yr'), N_t_ff_max*units['time_unit'].in_units('yr')]
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), total_stars)
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), N_vis_array)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("N$_{stars}$")
        plt.axhline(y=N_vis_val)
        plt.axvline(x=time_bounds[0])
        plt.axvline(x=time_bounds[1])
        plt.savefig("visible_stars.jpg")
elif args.match_method == 4:
    M_tot_value = 150
    M_tot = np.sum(global_data['m'], axis=1)*units['mass_unit'].value
    if args.use_t_spread == 'True':
        M_tot_ind = np.argmin(np.abs(M_tot_value-M_tot))
        M_tot_t_ff = global_data['time'].T[0][M_tot_ind]
        time_bounds = [(M_tot_t_ff-dt)*units['time_unit'].in_units('yr'),(M_tot_t_ff+dt)*units['time_unit'].in_units('yr')]
    else:
        M_tot_spread = M_tot_value*args.threshold_spread
        M_tot_min = np.argmin(np.abs((M_tot_value - M_tot_spread)-M_tot))
        M_t_ff_min = global_data['time'].T[0][M_tot_min]
        M_tot_max = np.argmin(np.abs((M_tot_value + M_tot_spread)-M_tot))
        M_t_ff_max = global_data['time'].T[0][M_tot_max]
        time_bounds = [M_t_ff_min*units['time_unit'].in_units('yr'), M_t_ff_max*units['time_unit'].in_units('yr')]
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), M_tot)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("M$_{total}$")
        plt.axhline(y=M_tot_value)
        plt.axvline(x=time_bounds[0])
        plt.axvline(x=time_bounds[1])
        plt.savefig("M_total.jpg")
elif args.match_method == 5:
    ratio_val = 1./3.
    R_prev = 0
    time_yrs = (global_data['time'].T[0] - global_data['time'].T[0][0])*units['time_unit'].in_units('yr')
    if args.use_t_spread == 'True':
        potential_inds = []
        Ratio_array = []
        for time in global_data['time'].T[0]:
            time_yt = yt.YTArray(time*scale_t, 's')
            global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
            sink_inds = np.where(global_data['m'][global_ind]>0)[0]
            L_tot = luminosity(global_data, sink_inds, global_ind)
            M_dot = accretion(sink_inds, global_ind)
            vis_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
            N_vis = len(vis_inds)
            ratio = N_vis/len(sink_inds)
            Ratio_array.append(ratio)
            if ratio == ratio_val:
                potential_inds.append(global_ind)
            elif ratio_val > R_prev and ratio_val < ratio:
                potential_inds.append(global_ind)
            R_prev = ratio
        if rank == 0:
            plt.clf()
            plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), Ratio_array)
            plt.xlabel("Time (t$_{ff}$)")
            plt.ylabel("Ratio $N_{vis}/N_{sinks}$")
            plt.axhline(y=ratio_val)
            plt.savefig("ratio.jpg")
        R_ind = int(len(potential_inds)/2)
        R_t_ff = global_data['time'].T[0][potential_inds[R_ind]]
        time_bounds = [(R_t_ff-dt)*units['time_unit'].in_units('yr'),(R_t_ff+dt)*units['time_unit'].in_units('yr')]
    else:
        ratio_spread = ratio_val*args.threshold_spread
        potential_min_inds = []
        potential_max_inds = []
        Ratio_array = []
        for time in global_data['time'].T[0]:
            time_yt = yt.YTArray(time*scale_t, 's')
            global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
            sink_inds = np.where(global_data['m'][global_ind]>0)[0]
            L_tot = luminosity(global_data, sink_inds, global_ind)
            M_dot = accretion(sink_inds, global_ind)
            vis_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
            N_vis = len(vis_inds)
            ratio = N_vis/len(sink_inds)
            Ratio_array.append(ratio)
            if ratio == (ratio_val-ratio_spread):
                potential_min_inds.append(global_ind)
            elif (ratio_val-ratio_spread)> R_prev and (ratio_val-ratio_spread)<ratio:
                potential_min_inds.append(global_ind)
            if ratio == (ratio_val+ratio_spread):
                potential_max_inds.append(global_ind)
            elif (ratio_val+ratio_spread) > R_prev and (ratio_val+ratio_spread)<ratio:
                potential_max_inds.append(global_ind)
            R_prev = ratio
        if rank == 0:
            plt.clf()
            plt.plot(global_data['time'].T[0].in_units('yr'), total_stars)
            plt.plot(global_data['time'].T[0].in_units('yr'), N_vis_array)
            plt.xlabel("Time (t$_{ff}$)")
            plt.ylabel("Ratio $N_{vis}/N_{sinks}$")
            plt.axhline(y=N_vis_val)
            plt.savefig("ratio.jpg")
        N_min = int(len(potential_min_inds)/2)
        N_max = int(len(potential_max_inds)/2)
        N_t_ff_min = global_data['time'].T[0][potential_min_inds[N_min]]
        try:
            N_t_ff_max = global_data['time'].T[0][potential_max_inds[N_max]]
        except:
            N_t_ff_max = global_data['time'].T[0][-1]
        time_bounds = [N_t_ff_min*units['time_unit'].in_units('yr'), N_t_ff_max*units['time_unit'].in_units('yr')]
    

try:
    start_time_ind = np.argmin(abs(global_data['time'].T[0]*units['time_unit'].in_units('yr')-time_bounds[0]))
    end_time_ind = np.argmin(abs(global_data['time'].T[0]*units['time_unit'].in_units('yr')-time_bounds[1]))
except:
    start_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[0]))
    end_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[1]))

if rank == 0:
    print("tstart, tend", time_bounds[0].in_units('kyr'), time_bounds[1].in_units('kyr'))
    print("SFE_start, SFE_end", (np.sum(global_data['m'][start_time_ind])*units['mass_unit'].value)/units['mass_unit'].value, (np.sum(global_data['m'][end_time_ind])*units['mass_unit'].value)/units['mass_unit'].value)
    print("nstars_start, nstars_end", len(np.where(global_data['m'][start_time_ind]>0)[0]), len(np.where(global_data['m'][end_time_ind]>0)[0]))
    print("mstars_start, mstars_end", np.sum(global_data['m'][start_time_ind])*units['mass_unit'].value, np.sum(global_data['m'][end_time_ind])*units['mass_unit'].value)
    print("number of records", end_time_ind-start_time_ind)

radius = yt.YTQuantity(2.0, 'rsun')
temperature = yt.YTQuantity(3000, 'K')

if args.starting_ind != None:
    update = True
    start_time_ind = args.starting_ind
elif args.make_plots_only != 'False':
    update = False
elif len(Times) == (end_time_ind-start_time_ind+1):
    update = False
else:
    update = True
    start_time_ind = start_time_ind + len(Times)
    
#============================================================================================
#Now we enter the actual multiplicity analysis
sys.stdout.flush()
CW.Barrier()

if args.debugging == "True":
    time_its = [10815, 10818, 10821, 10825]
else:
    time_its = range(start_time_ind, end_time_ind+1)

if update == True and args.make_plots_only == 'False':
    rit = -1
    for time_it in time_its:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            #"""
            n_stars = np.where(global_data['m'][time_it]>0)[0]
            abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][n_stars])
            time = global_data['time'][time_it][n_stars][0]
            sfe = np.sum(mass)
            
            sink_inds = np.where(global_data['m'][time_it]>0)[0]
            L_tot = luminosity(global_data, sink_inds, time_it)
            M_dot = accretion(sink_inds, time_it)
            vis_inds_tot = np.where((L_tot>=luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
            
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l.value
            S._scale_v = scale_v.value
            S._scale_t = scale_t.value
            S._scale_d = scale_d.value
            S._time = yt.YTArray(time, '')
            if args.visible_only == "True":
                S._abspos = yt.YTArray(abspos[vis_inds_tot], '')
                S._absvel = yt.YTArray(absvel[vis_inds_tot], '')
                S._mass = yt.YTArray(mass[vis_inds_tot], '')
            else:
                S._abspos = yt.YTArray(abspos, '')
                S._absvel = yt.YTArray(absvel, '')
                S._mass = yt.YTArray(mass, '')
            
            time_yt = yt.YTArray(time*scale_t, 's')
            Times.append(int(time_yt.in_units('yr').value))
            SFE.append(sfe)
            
            for sink_ind in sink_inds:
                if str(sink_ind) not in Sink_Luminosities.keys():
                    Sink_Luminosities.update({str(sink_ind): [[int(time_yt.in_units('yr').value), float(L_tot[sink_ind].value)]]})
                    Sink_Accretion.update({str(sink_ind): [[int(time_yt.in_units('yr').value), float(M_dot[sink_ind].value)]]})
                else:
                    Sink_Luminosities[str(sink_ind)].append([int(time_yt.in_units('yr').value), float(L_tot[sink_ind].value)])
                    Sink_Accretion[str(sink_ind)].append([int(time_yt.in_units('yr').value), float(M_dot[sink_ind].value)])
            
            CF_per_bin = []
            n_systems = []
            
            for bin_it in range(1,len(S_bins)):
                if len(n_stars) > 1:
                    if multiplicity_analysis_projection == False:
                        res = m.multipleAnalysis(S,cutoff=S_bins[bin_it], bound_check=bound_check, nmax=6, Grho=Grho, max_iter=100, use_mid_point_sep=use_mid_point_sep)#cyclic=False
                    else:
                        res = m.multipleAnalysis(S,cutoff=S_bins[bin_it], bound_check=bound_check, nmax=6, projection=multiplicity_analysis_projection, axis=args.axis, projection_vector=proj_unit, Grho=Grho, max_iter=100, use_mid_point_sep=use_mid_point_sep)#cyclic=False
                        if args.axis == 'x':
                            zero_ind = 0
                        elif args.axis == 'y':
                            zero_ind = 1
                        else:
                            zero_ind = 2
                        res['midpoint'].T[zero_ind] = 0
                        if len(np.where(res['midpointSep']>0)[0]) > 0:
                            update_midspoint_sep = np.where(res['midpointSep']>0)[0][0]
                            pos1 = res['midpoint'][res['index1'][update_midspoint_sep:]]
                            pos2 = res['midpoint'][res['index2'][update_midspoint_sep:]]
                            mid_sep = np.sqrt(np.sum(np.square(pos1 - pos2), axis=1))
                            res['midpointSep'][update_midspoint_sep:] = mid_sep
                    
                    #else:
                    #    res = m.multipleAnalysis(S,cutoff=S_bins[bin_it], bound_check=bound_check, nmax=6, projection=True, axis=args.axis, projection_vector=proj_unit, Grho=Grho)#cyclic=False
                else:
                      res = {'abspos'      : abspos*units['length_unit'].in_units('AU'),
                            'absvel'      : absvel,
                            'mass'        : mass*units['mass_unit'].in_units('Msun'),
                            'n'           : np.ones(mass.shape,dtype=np.int),
                            'ekin'        : np.ones(mass.shape,dtype=np.float),
                            'epot'        : np.ones(mass.shape,dtype=np.float),
                            'reducedMass' : np.ones(mass.shape,dtype=np.float),
                            'separation'  : np.zeros(mass.shape,dtype=np.float),
                            'separation_vector': np.zeros(abspos.shape,dtype=np.float),
                            'midpoint'    : abspos,
                            'midpointSep' : np.zeros(mass.shape,dtype=np.float),
                            'relativeSpeed' : np.zeros(mass.shape,dtype=np.float),
                            'semiMajorAxis' : np.zeros(mass.shape,dtype=np.float),
                            'eccentricity'  : np.ones(mass.shape,dtype=np.float),
                            'index1'      : np.zeros(mass.shape,dtype=np.int)-1,
                            'index2'      : np.zeros(mass.shape,dtype=np.int)-1,
                            'indices'     : np.arange(mass.size,dtype=np.int),
                            'topSystem'   : np.ones(mass.shape,dtype=np.bool),
                            'cutoff'      : S_bins[bin_it],
                            'nmax'        : 6,
                            'converged'   : False,
                            'time'        : time}
                
                sink_inds = np.where((res['n']==1))[0]
                sink_inds_total = np.arange(len(res['n']))
                nan_size = len(sink_inds_total) - len(sink_inds)
                nan_array = yt.YTArray(np.ones(nan_size)*np.nan, 'Lsun')
                #sink_inds = np.arange(len(res['n']))
                #If ONLY using visible stars, makes sure to recalculate the L_tot and M_tot for only those!
                if args.visible_only == "True":
                    L_tot = luminosity(global_data, vis_inds_tot, time_it)
                    M_dot = accretion(vis_inds_tot, time_it)
                    L_tot = np.append(L_tot, nan_array)
                    M_dot = np.append(M_dot, nan_array)
                    vis_inds = np.where((L_tot>=luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
                else:
                    L_tot = luminosity(global_data, sink_inds, time_it)
                    M_dot = accretion(sink_inds, time_it)
                    L_tot = np.append(L_tot, nan_array)
                    M_dot = np.append(M_dot, nan_array)
                    vis_inds = np.where((L_tot>=luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<=args.upper_L_limit))[0]
                visible_stars = sink_inds[vis_inds]
                visible_subcomps = visible_stars[np.where(res['topSystem'][visible_stars]==False)]
                checked_visible_inds = []
            
                top_inds = np.where(res['topSystem'])[0]
                
                if args.verbose_printing != 'False':
                    if size == 1:
                        print("=========================================================")
                        print("Doing time_it", time_it, "of", end_time_ind+1)
                        print("TRUE NUMBER OF VISIBLE STARS IS", str(len(visible_stars)))
                    else:
                        print_lines = ["========================================================="]
                        print_line = "Doing time_it " + str(time_it) + " of " + str(end_time_ind+1)
                        print_lines.append(print_line)
                        print_line = "TRUE NUMBER OF VISIBLE STARS IS " + str(len(visible_stars))
                        print_lines.append(print_line)
                else:
                    if bin_it == 1:
                        print("Doing time_it", time_it, "of", end_time_ind+1)

                #Find all singles and top systems with separations below the bin lower bound
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0] #These are true singles
                s_fake = np.where((res[sep_key]<S_bins[bin_it-1])&(res['topSystem']==True)&(res['n']!=1))[0] #These are Top systems whose largest separation is below the separatino bin. But these separations are calculated using the center of mass.

                if args.verbose_printing != 'False':
                    print_line = "AND", len(set(s_true).intersection(set(visible_stars))), "ARE VISIBLE SINGLE STARS"
                    if size == 1:
                        print(print_line)
                    else:
                        print_lines.append(print_line)
                        
                visible_singles = list(set(s_true).intersection(set(visible_stars)))
                checked_visible_inds = visible_singles
                invisible_singles = np.setdiff1d(s_true,visible_singles)
                if args.verbose_printing != 'False':
                    print_line = "UP TO " + str(int(S_bins[bin_it])) +"AU:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds]))
                    if size == 1:
                        print(print_line)
                    else:
                        print_lines.append(print_line)
                        
                res['n'][invisible_singles] = 0
                if args.verbose_printing != 'False':
                    print_line = "AFTER REMOVING " + str(len(invisible_singles)) + " INVISIBLE SINGLE STARS:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds]))
                    if size == 1:
                        print(print_line)
                    else:
                        print_lines.append(print_line)
                
                removed_sys = 0
                removed_stars = 0
                redefined_n = 0
                vis_stars_in_collapsed_systems = 0
                #Check fake single systems
                for fake in s_fake:
                    sys_comps = losi(fake, res)
                    sys_comps = sorted(flatten(sys_comps))
                    
                    vis_s_fake_inds = set(sys_comps).intersection(set(visible_stars))
                    
                    if len(vis_s_fake_inds) == 0:
                        res['n'][fake] = 0
                        removed_sys = removed_sys + 1
                        removed_stars = removed_stars + len(sys_comps)
                        L_tot[fake] = 0.0
                        M_dot[fake] = 0.0
                    else:
                        res['n'][fake] = 1
                        L_tot[fake] = np.sum(L_tot[list(vis_s_fake_inds)])
                        M_dot[fake] = np.sum(M_dot[list(vis_s_fake_inds)])
                        vis_stars_in_collapsed_systems = vis_stars_in_collapsed_systems + len(vis_s_fake_inds)
                        redefined_n = redefined_n + 1
                        checked_visible_inds = checked_visible_inds + list(vis_s_fake_inds)
                                
                if args.verbose_printing != 'False':
                    print_line = "AFTER REMOVING " + str(removed_sys) + " SYSTEMS WITH", str(removed_stars),"STARS, AND COLLAPSING " + str(redefined_n) + " SYSTEMS (WITH " + str(vis_stars_in_collapsed_systems) + " VISIBLE STARS):" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds]))
                    if size == 1:
                        print(print_line)
                    else:
                        print_lines.append(print_line)
                
                #Determine which systems could still be multiples
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True) & (res[sep_key]>S_bins[bin_it-1]))[0]
                
                #Let's go over the multiple systems and remove subsystems that are below the separation bin, or invisible
                removed_stars = 0
                removed_systems = [0, 0, 0, 0, 0, 0, 0]
                redefined_systems = [0, 0, 0, 0, 0, 0, 0]
                added_systems = [0, 0, 0, 0, 0, 0, 0]
                for multi_ind in multi_inds:
                    sys_comps = losi(multi_ind, res)
                    
                    vis_multi_ind = set(sorted(flatten(sys_comps))).intersection(set(visible_stars))
                    
                    #if len(vis_multi_ind) < len(flatten(sys_comps)):
                    redef_ind = len(flatten(sys_comps)) - 1
                    redefined_systems[redef_ind] = redefined_systems[redef_ind] + 1
                    if len(vis_multi_ind) == 0:
                        L_tot[multi_ind] = 0.0
                        M_dot[multi_ind] = 0.0
                        res['n'][multi_ind] = 0
                        removed_stars = removed_stars + len(flatten(sys_comps))
                    elif len(vis_multi_ind) == 1:
                        L_tot[multi_ind] = L_tot[list(vis_multi_ind)]
                        M_dot[multi_ind] = M_dot[list(vis_multi_ind)]
                        res['n'][multi_ind] = 1
                        removed_stars = removed_stars + len(flatten(sys_comps)) - 1
                        added_systems[0] = added_systems[0] + 1
                        checked_visible_inds = checked_visible_inds + list(vis_multi_ind)
                    elif len(vis_multi_ind) > 1:
                        checked_visible_inds = checked_visible_inds + list(vis_multi_ind)
                        invisible_components = list(set(flatten(sys_comps)).symmetric_difference(vis_multi_ind))
                        #Remove invisible inds
                        system_structure = losi(multi_ind, res)
                        sys_string = str(system_structure)
                        
                        reduced = False
                        while reduced == False:
                            bracket_pos = []
                            for char_it in range(len(sys_string)):
                                if sys_string[char_it] == '[':
                                    bracket_pos.append(char_it)
                                elif sys_string[char_it] == ']':
                                    open_ind = bracket_pos.pop()
                                    try:
                                        sub_sys_comps = eval(sys_string[open_ind:char_it+1])
                                    except:
                                        comma_split = sys_string[open_ind:char_it+1].split(',')
                                        sub_sys_comps = eval(comma_split[0]+comma_split[1])
                                    if len(sub_sys_comps) == 2:
                                        binary_ind = np.where((res['index1']==sub_sys_comps[0])&(res['index2']==sub_sys_comps[1]))[0][0]
                                        ind_1 = res['index1'][binary_ind]
                                        ind_2 = res['index2'][binary_ind]
                                        if use_mid_point_sep:
                                            pos_diff = res['midpoint'][ind_1] - res['midpoint'][ind_2]
                                        else:
                                            pos_diff = res['separation_vector'][ind_1] - res['separation_vector'][ind_2]
                                        #MAKE SURE 2D SEP IS BEING SAVED.
                                        sep_value = np.sqrt(np.sum(pos_diff**2))
                                        if sep_value > (scale_l.in_units('AU')/2):
                                            update_inds = np.where(abs(pos_diff)>scale_l.in_units('AU')/2)[0]
                                            for ind in update_inds:
                                                if pos_diff[ind] < 0:
                                                    pos_diff[ind] = pos_diff[ind] + scale_l.in_units('AU').value
                                                else:
                                                    pos_diff[ind] = pos_diff[ind] - scale_l.in_units('AU').value
                                            sep_value = np.sqrt(np.sum(pos_diff**2))
                                            if sep_value > (scale_l.in_units('AU')/2):
                                                print("FAILED ON time_it:", time_it, "SEPARATION=", sep_value)
                                                if args.debugging == "True":
                                                    import pdb
                                                    pdb.set_trace()
                                        if sep_value < S_bins[bin_it-1]:
                                            #Reduce systems! Let's check if any of the components are visible
                                            vis_subs = set([ind_1, ind_2]).intersection(set(visible_stars))
                                            if len(vis_subs) > 0:
                                                L_tot[binary_ind] = np.max(L_tot[list(vis_subs)]) #np.sum(L_tot[list(vis_subs)])
                                                M_dot[binary_ind] = np.max(M_dot[list(vis_subs)]) #np.sum(M_dot[list(vis_subs)])
                                                if use_mid_point_sep:
                                                    res['midpoint'][binary_ind] = (res['midpoint'][ind_1] + res['midpoint'][ind_2])/2#res['midpoint'][central_ind]
                                                else:
                                                    res['separation_vector'][binary_ind] = (res['separation_vector'][ind_1]*res['mass'][ind_1] + res['separation_vector'][ind_2]*res['mass'][ind_2])/(res['mass'][ind_1]+res['mass'][ind_2])
                                                replace_string = str(binary_ind)
                                                res['n'][multi_ind] = res['n'][multi_ind] - 1
                                                removed_stars = removed_stars + 1
                                            else:
                                                L_tot[binary_ind] = 0.0
                                                M_dot[binary_ind] = 0.0
                                                replace_string = ""
                                                res['n'][multi_ind] = res['n'][multi_ind] - 2
                                                removed_stars = removed_stars + 2
                                        else:
                                            vis_subs = set([ind_1, ind_2]).intersection(set(visible_stars))
                                            if len(vis_subs) > 0:
                                                L_tot[binary_ind] = np.max(L_tot[list(vis_subs)])#np.sum(L_tot[list(vis_subs)])
                                                M_dot[binary_ind] = np.max(M_dot[list(vis_subs)])#np.sum(M_dot[list(vis_subs)])
                                                if use_mid_point_sep:
                                                    res['midpoint'][binary_ind] = (res['midpoint'][ind_1] + res['midpoint'][ind_2])/2 #res['midpoint'][central_ind]
                                                else:
                                                    res['separation_vector'][binary_ind] = (res['separation_vector'][ind_1]*res['mass'][ind_1] + res['separation_vector'][ind_2]*res['mass'][ind_2])/(res['mass'][ind_1]+res['mass'][ind_2])
                                                replace_string = str(binary_ind)
                                                if len(vis_subs) == 1:
                                                    res['n'][multi_ind] = res['n'][multi_ind] - 1
                                                    removed_stars = removed_stars + 1
                                            else:
                                                L_tot[binary_ind] = 0.0
                                                M_dot[binary_ind] = 0.0
                                                replace_string = ""
                                                res['n'][multi_ind] = res['n'][multi_ind] - 2
                                                removed_stars = removed_stars + 2
                                        str_1 = sys_string[:open_ind]
                                        str_2 = sys_string[char_it+1:]
                                        sys_string = str_1 + replace_string + str_2
                                        if '[' not in sys_string:
                                            reduced = True
                                        break
                                    elif len(sub_sys_comps) == 1:
                                        binary_ind = np.where((res['index1']==sub_sys_comps[0])|(res['index2']==sub_sys_comps[0]))[0][0]
                                        vis_subs = set(sub_sys_comps).intersection(set(visible_stars))
                                        L_tot[binary_ind] = np.sum(L_tot[list(vis_subs)])
                                        M_dot[binary_ind] = np.sum(M_dot[list(vis_subs)])
                                        if len(vis_subs)>0:
                                            if use_mid_point_sep:
                                                res['midpoint'][binary_ind] = res['abspos'][list(vis_subs)][0]
                                            else:
                                                res['separation_vector'][binary_ind] = res['abspos'][list(vis_subs)][0]
                                            replace_string = str(binary_ind)
                                        else:
                                            replace_string = ""
                                            res['n'][multi_ind] = res['n'][multi_ind] - 1
                                            removed_stars = removed_stars + 1
                                        str_1 = sys_string[:open_ind]
                                        str_2 = sys_string[char_it+1:]
                                        sys_string = str_1 + replace_string + str_2
                                        if '[' not in sys_string:
                                            reduced = True
                                        break
                                    else:
                                        replace_string = ""
                                        str_1 = sys_string[:open_ind]
                                        str_2 = sys_string[char_it+1:]
                                        sys_string = str_1 + replace_string + str_2
                                        if '[' not in sys_string:
                                            reduced = True
                                        break
                        add_ind = res['n'][multi_ind] - 1
                        added_systems[add_ind] = added_systems[add_ind] + 1
                
                if args.verbose_printing != 'False':
                    print_line = "AFTER REDEFINING", str(redefined_systems), "SYSTEMS WITH", str(removed_stars),"INVISIBLE COMPONENTS TO", str(added_systems), ":" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds]))
                    if size == 1:
                        print(print_line)
                        print("TOTAL NUMBER OF STARS =", str(np.sum(res['n'][top_inds])))
                    else:
                        print_lines.append(print_line)
                        print_line = "TOTAL NUMBER OF STARS = " + str(np.sum(res['n'][top_inds]))
                        print_lines.append(print_line)
                    
                ns = len(np.where(res['n'][top_inds]==1)[0])
                nb = len(np.where(res['n'][top_inds]==2)[0])
                nt = len(np.where(res['n'][top_inds]==3)[0])
                nq = len(np.where(res['n'][top_inds]==4)[0])
                nq5 = len(np.where(res['n'][top_inds]==5)[0])
                ns6 = len(np.where(res['n'][top_inds]==6)[0])
                ns7 = len(np.where(res['n'][top_inds]==7)[0])
                n_systems.append([ns,nb,nt,nq,nq5,ns6,ns7])
                if (ns+nb+nt+nq+nq5+ns6+ns7) == 0:
                    cf = 0
                else:
                    cf = np.sum(np.arange(7)*np.array(n_systems[-1]))/np.sum(n_systems[-1])
                if args.verbose_printing != 'False':
                    print_line = "CF = " + str(cf)
                    if size == 1:
                        print("CF =", cf)
                    else:
                        print_lines.append(print_line)
                CF_per_bin.append(cf)
                '''
                if bin_it == len(S_bins)-1 and cf>0.5:
                    import pdb
                    pdb.set_trace()
                '''
                
                if len(vis_inds) == len(checked_visible_inds):
                    if args.verbose_printing != 'False':
                        if size == 1:
                            print("All visible stars are accounted for")
                        else:
                            print_lines.append("All visible stars are accounted for")
                else:
                    print("NOT ALL VISIBLE STARS HAVE BEEN ACCOUNTED FOR!")
                    import pdb
                    pdb.set_trace()
                if size > 1 and args.verbose_printing != 'False':
                    for print_line in print_lines:
                        print(print_line)

            CF_arrays.append(CF_per_bin)
            N_sys_total.append(n_systems)
            pickle_file_rank = pickle_file.split('.pkl')[0] + "_" +str(rank) + ".pkl"
            file = open(pickle_file_rank, 'wb')
            pickle.dump((CF_arrays, N_sys_total, Times, SFE, Sink_Luminosities, Sink_Accretion),file)
            file.close()
            print('updated pickle', pickle_file_rank, "for time_it", time_it, "of", end_time_ind+1)
                
    sys.stdout.flush()
    CW.Barrier()
    print("FINISHED GOING THROUGH TIMES ON RANK", rank)
    
    
sys.stdout.flush()
CW.Barrier()
#=====================================================
#Create plots below

if rank == 0:
    #compile together data
    try:
        file = open(pickle_file+'.pkl', 'rb')
        Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion = pickle.load(file)
        file.close()
    except:
        pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
        Times_full = []
        SFE_full = []
        CF_per_bin_full = []
        n_systems_full = []
        for pick_file in pickle_files:
            file = open(pick_file, 'rb')
            CF_arrays, N_sys_total, Times, SFE, Sink_Luminosities, Sink_Accretion = pickle.load(file)
            file.close()
            Sink_Luminosities_full = {}
            Sink_Accretion_full = {}
            for L_key in Sink_Luminosities.keys():
                if L_key not in Sink_Luminosities_full.keys():
                    Sink_Luminosities_full.update({L_key:Sink_Luminosities[L_key]})
                    Sink_Accretion_full.update({L_key:Sink_Accretion[L_key]})
                else:
                    Sink_Luminosities_full[L_key] = Sink_Luminosities_full[L_key] + Sink_Luminosities[L_key]
                    Sink_Accretion_full[L_key] = Sink_Accretion_full[L_key] + Sink_Accretion[L_key]
            Times_full = Times_full + Times
            SFE_full = SFE_full + SFE
            CF_per_bin_full = CF_per_bin_full + CF_arrays
            n_systems_full = n_systems_full + N_sys_total
            os.remove(pick_file)
        
        #Let's sort the data
        for L_key in Sink_Luminosities_full:
            sorted_array = np.array(Sink_Luminosities_full[L_key])
            dict_sort = np.argsort(sorted_array.T[0])
            sorted_array.T[0] = np.array(Sink_Luminosities_full[L_key]).T[0][dict_sort]
            sorted_array.T[1] = np.array(Sink_Luminosities_full[L_key]).T[1][dict_sort]
            Sink_Luminosities[L_key] = sorted_array.tolist()
            sorted_array.T[0] = np.array(Sink_Accretion_full[L_key]).T[0][dict_sort]
            sorted_array.T[1] = np.array(Sink_Accretion_full[L_key]).T[1][dict_sort]
            Sink_Accretion[L_key] = sorted_array.tolist()
        sorted_inds = np.argsort(Times_full)
        Times = np.array(Times_full)[sorted_inds]
        SFE = np.array(SFE_full)[sorted_inds]
        CF_Array_Full = np.array(CF_per_bin_full)[sorted_inds]
        N_sys_total = np.array(n_systems_full)[sorted_inds]
        
        file = open(pickle_file+'.pkl', 'wb')
        pickle.dump((Times, SFE, CF_Array_Full, N_sys_total, Sink_Luminosities, Sink_Accretion),file)
        file.close()
    
    Summed_systems = np.sum(np.array(N_sys_total), axis=0)
    CF_top = Summed_systems[:,1] + Summed_systems[:,2]*2 + Summed_systems[:,3]*3 + Summed_systems[:,4]*4 + Summed_systems[:,5]*5 + Summed_systems[:,6]*6
    CF_bot = np.sum(Summed_systems, axis=1)
    CF_Total = CF_top/CF_bot

    file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data.pkl", "rb")
    bin_centers, CF_per_bin_Tobin_Per = pickle.load(file_open)
    file_open.close()

    file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Orion_data.pkl", "rb")
    bin_centers, CF_per_bin_Tobin_Ori = pickle.load(file_open)
    file_open.close()


    plt.clf()
    try:
        plt.bar(bin_centers, CF_Total, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
        plt.bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', alpha=0.5, label="Perseus")
        plt.bar(bin_centers, CF_per_bin_Tobin_Ori, width=0.25, edgecolor='black', alpha=0.5, label="Orion")
    except:
        plt.bar(bin_centers[1:], CF_Total, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
        plt.bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', alpha=0.5, label="Perseus")
        plt.bar(bin_centers, CF_per_bin_Tobin_Ori, width=0.25, edgecolor='black', alpha=0.5, label="Orion")
        S_bins = S_bins[1:]
    plt.legend(loc='best')
    plt.xlabel('Log Separation (AU)')
    plt.ylabel('Companion Frequency')
    plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
    plt.ylim(bottom=0.0)
    plt.savefig(savedir +  args.figure_prefix + 'Total_companion_frequency_.jpg', format='jpg', bbox_inches='tight')
    print('created Total_companion_frequency.jpg')
    
    #Create mean CF
    CF_median = []
    CF_err = []
    N = np.shape(N_sys_total)[0]

    #plt.clf()
    #fig, axs = plt.subplots(3, 4, constrained_layout=True, sharex=True, sharey=True) #figsize=(4, 3*len(pickle_files))
    #axs = axs.flatten()
    #plt.subplots_adjust(wspace=0.0, hspace=0.0)
    
    plt.clf()
    fig = plt.figure()
    total_number_of_figures = len(bin_centers)
    columns = 4
    rows = int(np.ceil(total_number_of_figures/columns))
    fig.set_size_inches(3.25*columns, 3.25*rows)
    gs = gridspec.GridSpec(rows, columns)\
    
    gs.update(wspace=0.0, hspace=0.0)
    axes_dict = {}

    for bit in range(len(S_bins)-1):
        ax_label = 'ax' + str(bit)
        if bit == 0:
            axes_dict.update({ax_label:fig.add_subplot(gs[0,bit])})
            xticklabels = axes_dict[ax_label].get_xticklabels()
            plt.setp(xticklabels, visible=False)
            axes_dict[ax_label].tick_params(axis="x",direction="in")
            axes_dict[ax_label].set_xlim([0, 0.4])
            axes_dict[ax_label].set_ylim([0, 1600])
            #axes_dict[ax_label].set_ylim(bottom=0.0)
            #axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
            axes_dict[ax_label].set_ylabel("Count")
            axes_dict[ax_label].set_xlabel("CF")
            axes_dict[ax_label].tick_params(axis="y",direction="in")
        else:
            axes_dict.update({ax_label:fig.add_subplot(gs[int(bit/columns),np.remainder(bit,columns)], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
            if bit > 6:
                xticklabels = axes_dict[ax_label].get_xticklabels()
                plt.setp(xticklabels[0], visible=False)
                if bit == 8:
                    plt.setp(xticklabels[0], visible=True)
            else:
                xticklabels = axes_dict[ax_label].get_xticklabels()
                plt.setp(xticklabels, visible=False)
            if np.remainder(bit,columns) != 0:
                yticklabels = axes_dict[ax_label].get_yticklabels()
                plt.setp(yticklabels, visible=False)
                yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
                plt.setp(yticklabels, visible=False)
            else:
                #axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                axes_dict[ax_label].set_ylabel("Count")
            axes_dict[ax_label].set_xlabel("CF")
            axes_dict[ax_label].tick_params(axis="x",direction="in")
            axes_dict[ax_label].tick_params(axis="y",direction="in")
        time_text = plt.text(0.08, 1500, 'Log Separtion (AU) = ['+str(np.log10(S_bins)[bit])+','+str(np.log10(S_bins)[bit+1])+']', va="center", ha="left", color='k')
        
        #err_bins = np.linspace(0, 0.4, 41)
        err_bins = np.linspace(0, 0.4, 21)
        #err_bins = np.linspace(np.min(np.array(CF_Array_Full)[:,bit]), np.max(np.array(CF_Array_Full)[:,bit]), 15)
        err_centers = (err_bins[:-1]+err_bins[1:])/2
        err_hist, err_bins = np.histogram(CF_Array_Full[:,bit], bins=err_bins)
        width = err_bins[1] - err_bins[0]
        axes_dict[ax_label].bar(err_centers, err_hist, width=width)
        #plt.bar(err_centers, err_hist, width=width)
        
        xs = np.ones(np.shape(CF_Array_Full[:,bit]))*bin_centers[bit]
        #plt.scatter(xs, np.array(CF_Array_Full)[:,bit], alpha=0.5, color='b')
        median = np.median(CF_Array_Full[:,bit])
        axes_dict[ax_label].axvline(x=median, label='median', color='r', alpha=0.25)
        #plt.axvline(x=median, label='median', color='r')
        #plt.scatter(bin_centers[bit], median, color='r')
        mean = np.mean(CF_Array_Full[:,bit])
        axes_dict[ax_label].axvline(x=mean, label='mean', color='orange', alpha=0.25)
        #plt.axvline(x=mean, label='mean', color='orange')
        #plt.scatter(bin_centers[bit], median, color='orange')
        std = np.std(CF_Array_Full[:,bit], ddof=1)
        standard_error = std/np.sqrt(N)
        #rough_and_ready_err = (np.max(np.array(CF_Array_Full)[:,bit]) - np.min(np.array(CF_Array_Full)[:,bit]))*(2./3.)
        axes_dict[ax_label].axvline(x=mean-std, label='std', color='b', alpha=0.25)
        axes_dict[ax_label].axvline(x=mean+std, label='std', color='b', alpha=0.25)
        #plt.axvline(x=mean-std, label='std', color='b')
        #plt.axvline(x=mean+std, label='std', color='b')
        #plt.axvline(x=mean-standard_error, label='standard error', color='b', linestyle=":")
        #plt.axvline(x=mean+standard_error, label='standard error', color='b', linestyle=":")
        if bit == 10:
            axes_dict[ax_label].legend(loc='best')#,bbox_to_anchor=(0.985, 0.5))
        #plt.ylim([0, 1600])
        #plt.savefig(savedir+"error_dist_bin_"+str(bit)+".jpg", format='jpg', bbox_inches='tight')
        standard_deviation = [median-(mean-std), (mean+std)-median]
        #plt.plot([bin_centers[bit], bin_centers[bit]], [mean-std, mean+std], color='black')
        CF_median.append(median)
        CF_err.append(standard_deviation)

    
    plt.savefig(savedir+"error_dist.pdf", format='pdf', bbox_inches='tight')
    #plt.ylim(bottom=0.0)
    #plt.savefig(savedir+"error_dist.jpg")
    CF_err = np.array(CF_err)

    plt.clf()
    try:
        plt.bar(bin_centers, CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
    except:
        plt.bar(bin_centers[1:], CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
    plt.bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', alpha=0.5, label="Perseus")
    plt.bar(bin_centers, CF_per_bin_Tobin_Ori, width=0.25, edgecolor='black', alpha=0.5, label="Orion")
    plt.legend(loc='best')
    plt.xlabel('Log Separation (AU)')
    plt.ylabel('Companion Frequency')
    plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
    plt.ylim(bottom=0.0)
    plt.title("Median over " + str(np.shape(CF_Array_Full)[0]) + " histograms")
    plt.savefig(savedir+'Median_companion_frequency.pdf', format='pdf', bbox_inches='tight')
    print('created Median_companion_frequency.pdf')
    """
    #CF calculated by SUMMING all CF historgrams
    CF_Array_Full = np.array(CF_Array_Full)
    CF_sum = np.sum(CF_Array_Full, axis=0)

    plt.clf()
    plt.bar(bin_centers, CF_sum, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
    plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
    plt.legend(loc='best')
    plt.xlabel('Log Separation (AU)')
    plt.ylabel('Companion Frequency')
    plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
    plt.ylim(bottom=0.0)
    plt.savefig(savedir + args.figure_prefix + 'Sum_companion_frequency_.jpg', format='jpg', bbox_inches='tight')
    print('created Sum_companion_frequency.jpg')
    
    Sink_L_means = []
    for L_key in Sink_Luminosities.keys():
        Sink_L_means.append(np.mean(np.array(Sink_Luminosities[L_key]).T[1]))
    L_bins = np.logspace(-2,4,25)
    L_hist, L_bins = np.histogram(Sink_L_means, bins=L_bins)
    
    Tobin_luminosities_0_1 = np.array([1.8, 0.50, 0.6, 0.7, 0.4, 0.36, 1.40, 0.80, 0.43, 1.20, 3.70, 1.70, 0.16, 0.05, 1.6, 0.54, 0.3, 1.20, 23.2, 0.16, 4.7, 16.80, 0.54, 0.09, 0.24, 1.80, 1.90, 3.20, 0.69, 1.50, 0.90, 1.30, 4.20, 3.60, 8.40, 0.68, 2.6, 1.80, 0.40, 0.70, 0.30, 0.60, 19.00, 5.30, 1.50, 0.30, 0.50, 0.54, 0.17, 0.32, 0.70, 2.80, 6.9, 1.10, 4.00, 7.00, 0.10, 32.50, 1.00, 8.3, 9.2, 1.4, 0.87, 1.50, 3.20, 9.1, 0.63, 0.16])
    Tobin_hist, L_bins = np.histogram(Tobin_luminosities_0_1, bins=L_bins)
    
    plt.clf()
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Simulation")
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Tobin et al (2016)")
    plt.legend(loc='best')
    plt.xlabel('Mean Luminosty (log(L))')
    plt.ylabel('Number')
    plt.ylim(bottom=0.0)
    plt.savefig(savedir + args.figure_prefix + 'L_mean_hist.jpg', format='jpg', bbox_inches='tight')
    print('created L_mean_hist.jpg')
    """
    
print("FINISHED ON RANK", rank)
sys.stdout.flush()
CW.Barrier()
