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
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.09)
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickles", type=str, default='')
    parser.add_argument("-use_midpoint", "--use_midpoint_separation", help="Do you want to use the midpoint separation instread of separation", type=str, default='')
    parser.add_argument("-lifetime", "--lifetime_threshold", help="What life time threshold do you want to consider when making Luminosity histogram", type=float, default=10000)
    parser.add_argument("-proj_vec", "--projection_vector", help="What projection vector do you want to use?", type=str, default='')
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=10000)
    parser.add_argument("-use_t_s", "--use_t_spread", help="Do you want to define the time spread using t_spread instead of a spread defined by the method?", type=str, default='False')
    parser.add_argument("-match_meth", "--match_method", help="How do you want to select times? 1. SFE, 2.SPE/t_ff, 3. No. visible stars, 4. Total Accreted Mass", type=int, default=1)
    parser.add_argument("-start_ind", "--starting_ind", help="Do you want to start the analysis at a particular starting ind?", type=int, default=None)
    parser.add_argument("-n_vis_thres", "--n_visible_threshold", help="what threshold do you want to use for number fo stars?", type=int, default=106)
    parser.add_argument("-SFE_thres", "--SFE_threshold", help="What threshold do you want to use for the SFE?", type=float, default=0.04)
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

#What do I need to save?
Times = []
SFE = []
N_sys_total = []
CF_arrays = []
Sink_Luminosities = {}
Sink_Accretion = {}
All_separations = []

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
    SFE = np.sum(global_data['m'], axis=1)
    if args.use_t_spread == 'True':
        SFE_ind = np.argmin(np.abs(SFE_value-SFE))
        SFE_t_ff = global_data['time'].T[0][SFE_ind]*units['time_unit'].in_units('yr')
        time_bounds = [(SFE_t_ff-dt),(SFE_t_ff+dt)]
    else:
        SFE_spread = args.threshold_spread
        SFE_min_ind = np.argmin(np.abs((SFE_value-SFE_spread)-SFE))
        SFE_t_ff_min = global_data['time'].T[0][SFE_min_ind]*units['time_unit'].in_units('yr')
        SFE_max_ind = np.argmin(np.abs((SFE_value+SFE_spread)-SFE))
        SFE_t_ff_max = global_data['time'].T[0][SFE_max_ind]*units['time_unit'].in_units('yr')
        time_bounds = [SFE_t_ff_min,SFE_t_ff_max]
    if rank == 0:
        plt.clf()
        plt.plot(global_data['time'].T[0]*units['time_unit'].in_units('yr'), SFE)
        plt.xlabel("Time (t$_{ff}$)")
        plt.ylabel("SFE")
        plt.axhline(y=SFE_value)
        plt.axvline(x=time_bounds[0])
        plt.axvline(x=time_bounds[1])
        plt.savefig("SFE.jpg")
    del SFE
    SFE = []
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
    del SFE
    SFE = []
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

if args.update_pickle == 'True':
    update = True
elif args.update_pickle == 'False':
    update = False
    
#============================================================================================
#Now we enter the actual multiplicity analysis
sys.stdout.flush()
CW.Barrier()

time_its = range(start_time_ind, end_time_ind+1)
All_YSO_dens = []

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
            YSO_densities = []
            
            #For each visible star, find the distance to the 11th neighbour, then calculate the circular area (4pi*r^2), and the dnesity is 10/area.
            for vis_ind in vis_inds_tot:
                dx = abspos[vis_ind][0] - abspos.T[0]
                dy = abspos[vis_ind][1] - abspos.T[1]
                dz = abspos[vis_ind][2] - abspos.T[2]
                separation = np.sqrt(dx**2 + dy**2 + dz**2)
                neighbour_11 = np.sort(separation)[11]*units['length_unit'].in_units('AU')
                area = 4*np.pi*(neighbour_11**2)
                yso_dens = 10/area
                
                YSO_densities.append(yso_dens)
            Times.append(time)
            All_YSO_dens.append(YSO_densities)

            file = open('yso_dens_'+str(rank)+'.pkl', 'wb')
            pickle.dump((Time, All_YSO_dens), file)
            file.close()
            print("calculated YSO dens for time_it", time_it, "of", time_its[-1], ". Wrote file", 'yso_dens_'+str(rank)+'.pkl')

if rank == 0:
    pickle_files = sorted(glob.glob('yso_dens_*.pkl'))
    Times_full = []
    All_YSO_dens_full = []
    for pick_file in pickle_files:
        file = open(pick_file, 'rb')
        Time, All_YSO_dens = pickle.load(file)
        file.close()
        Times_full = Times_full + Times
        All_YSO_dens_full = All_YSO_dens_full + All_YSO_dens
        os.remove(pick_file)
            
    sorted_inds = np.argsort(Times_full)
    Times = np.array(Times_full)[sorted_inds]
    All_YSO_dens = np.array(All_YSO_dens_full)[sorted_inds]
    
    file = open('yso_dens_.pkl', 'wb')
    pickle.dump((Time, All_YSO_dens), file)
    file.close()
    
#Make YSO CDFs

