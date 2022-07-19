import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import pickle
import os
from mpi4py.MPI import COMM_WORLD as CW

f_acc= 0.5

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
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

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()

datadir = sys.argv[1]
savedir = sys.argv[2]
pickle_file = savedir + args.pickled_file
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

sys.stdout.flush()
CW.Barrier()

#What do I need to save?
Times = []
SFE = []
MF_true = []
MF_acc_lim = []
MF_acc_L_lower = []
MF_acc_L_upper = []
Single_L_mean = []
Single_L_std = []
Single_M_dot_mean = []
Single_M_dot_std = []

file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
dm = global_data['dm']*units['mass_unit'].in_units('Msun')
dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')
Accretion_array = dm/dt
print('Loaded global pickle data')

#dt == integration window, if you don't want to integrate over the entire simulation
time_bounds = [global_data['time'].T[0][0]*units['time_unit'].in_units('yr'), global_data['time'].T[0][-1]*units['time_unit'].in_units('yr')]
try:
    start_time_ind = np.argmin(abs(global_data['time'].T[0]*units['time_unit'].in_units('yr')-time_bounds[0]))
    end_time_ind = np.argmin(abs(global_data['time'].T[0]*units['time_unit'].in_units('yr')-time_bounds[1]))
except:
    start_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[0]))
    end_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[1]))

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
            if len(n_stars) > 1:
                abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                mass = np.array(global_data['m'][time_it][n_stars])
                time = global_data['time'][time_it][n_stars][0]
                sfe = np.sum(mass)
                
                time_yt = yt.YTArray(time*scale_t, 's')
                Times.append(int(time_yt.in_units('yr').value))
                SFE.append(sfe)
                
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
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                
                L_tot = luminosity(global_data, n_stars, time_it)
                M_dot = accretion(n_stars, time_it)
                Single_L_mean.append(np.mean(L_tot[s_true]))
                Single_L_std.append(np.std(L_tot[s_true]))
                Single_M_dot_mean.append(np.mean(M_dot[s_true]))
                Single_M_dot_std.append(np.std(M_dot[s_true]))
                
                total_systems = len(s_true) + len(multi_inds)
                MF_value = len(multi_inds)/total_systems
                MF_true.append(MF_value)
                
                #Multiplicity with accretion limit
                vis_inds_tot = np.where((M_dot>args.accretion_limit))[0]
                
                S._abspos = yt.YTArray(abspos[vis_inds_tot], '')
                S._absvel = yt.YTArray(absvel[vis_inds_tot], '')
                S._mass = yt.YTArray(mass[vis_inds_tot], '')
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, Grho=Grho, max_iter=100)
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                
                total_systems = len(s_true) + len(multi_inds)
                MF_value = len(multi_inds)/total_systems
                MF_acc_lim.append(MF_value)
                
                #Multiplicity with accretion and lower luminosity limit
                vis_inds_tot = np.where((L_tot>=args.lower_L_limit)&(M_dot>args.accretion_limit))[0]
                
                S._abspos = yt.YTArray(abspos[vis_inds_tot], '')
                S._absvel = yt.YTArray(absvel[vis_inds_tot], '')
                S._mass = yt.YTArray(mass[vis_inds_tot], '')
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, Grho=Grho, max_iter=100)
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                
                total_systems = len(s_true) + len(multi_inds)
                MF_value = len(multi_inds)/total_systems
                MF_acc_L_lower.append(MF_value)
                
                #Multiplicity with accretion and lower and upper luminosity limit
                vis_inds_tot = np.where((L_tot>=args.lower_L_limit)&(M_dot>args.accretion_limit)&(L_tot<=args.upper_L_limit))[0]
                
                S._abspos = yt.YTArray(abspos[vis_inds_tot], '')
                S._absvel = yt.YTArray(absvel[vis_inds_tot], '')
                S._mass = yt.YTArray(mass[vis_inds_tot], '')
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, Grho=Grho, max_iter=100)
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
                multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
                
                total_systems = len(s_true) + len(multi_inds)
                MF_value = len(multi_inds)/total_systems
                MF_acc_L_upper.append(MF_value)
                
                pickle_file_rank = pickle_file.split('.pkl')[0] + "_" +str(rank) + ".pkl"
                file = open(pickle_file_rank, 'wb')
                pickle.dump((Times, SFE, MF_true, MF_acc_lim, MF_acc_L_lower, MF_acc_L_upper, Single_L_mean, Single_L_std, Single_M_dot_mean, Single_M_dot_std),file)
                file.close()
                print('updated pickle', pickle_file_rank, "for time_it", time_it, "of", end_time_ind+1)
            
sys.stdout.flush()
CW.Barrier()
print("FINISHED GOING THROUGH TIMES ON RANK", rank)
    
if rank == 0:
    #compile together data
    try:
        file = open(pickle_file+'.pkl', 'rb')
        Times, SFE, MF_true, MF_acc_lim, MF_acc_L_lower, MF_acc_L_upper, Single_L_mean, Single_L_std, Single_M_dot_mean, Single_M_dot_std = pickle.load(file)
        file.close()
    except:
        pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
        Times_full = []
        SFE_full = []
        MF_true_full = []
        MF_acc_lim_full = []
        MF_acc_L_lower_full = []
        MF_acc_L_upper_full = []
        Single_L_mean_full = []
        Single_L_std_full = []
        Single_M_dot_mean_full = []
        Single_M_dot_std_full = []
        for pick_file in pickle_files:
            file = open(pick_file, 'rb')
            Times, SFE, MF_true, MF_acc_lim, MF_acc_L_lower, MF_acc_L_upper, Single_L_mean, Single_L_std, Single_M_dot_mean, Single_M_dot_std = pickle.load(file)
            file.close()
            Times_full = Times_full + Times
            SFE_full = SFE_full + SFE
            MF_true_full = MF_true_full + MF_true
            MF_acc_lim_full = MF_acc_lim_full + MF_acc_lim
            MF_acc_L_lower_full = MF_acc_L_lower_full + MF_acc_L_lower
            MF_acc_L_upper_full = MF_acc_L_upper_full + MF_acc_L_upper
            Single_L_mean_full = Single_L_mean_full + Single_L_mean
            Single_L_std_full = Single_L_std_full + Single_L_std
            Single_M_dot_mean_full = Single_M_dot_mean_full + Single_M_dot_mean
            Single_M_dot_std_full = Single_M_dot_std_full + Single_M_dot_std
            os.remove(pick_file)
        
        #Let's sort the data
        sorted_inds = np.argsort(Times_full)
        Times = np.array(Times_full)[sorted_inds]
        SFE = np.array(SFE_full)[sorted_inds]
        MF_true = np.array(MF_true_full)[sorted_inds]
        MF_acc_lim = np.array(MF_acc_lim_full)[sorted_inds]
        MF_acc_L_lower = np.array(MF_acc_L_lower_full)[sorted_inds]
        MF_acc_L_upper = np.array(MF_acc_L_upper_full)[sorted_inds]
        Single_L_mean = np.array(Single_L_mean_full)[sorted_inds]
        Single_L_std = np.array(Single_L_std_full)[sorted_inds]
        Single_M_dot_mean = np.array(Single_M_dot_mean_full)[sorted_inds]
        Single_M_dot_std = np.array(Single_M_dot_std_full)[sorted_inds]
        
        file = open(pickle_file+'.pkl', 'wb')
        pickle.dump((Times, SFE, MF_true, MF_acc_lim, MF_acc_L_lower, MF_acc_L_upper, Single_L_mean, Single_L_std, Single_M_dot_mean, Single_M_dot_std),file)
        file.close()
        
        plt.clf()
        plt.plot(SFE, MF_true, alpha=0.5, label="True")
        plt.plot(SFE, MF_acc_lim, alpha=0.5, label="Accretion limit")
        plt.plot(SFE, MF_acc_L_lower, alpha=0.5, label="Accretion + lower Luminosity limit")
        plt.plot(SFE, MF_acc_L_upper, alpha=0.5, label="Accretion + lower + upper Luminosity limit")
        plt.legend(loc='best')
        plt.xlabel('SFE')
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.ylabel('MF')
        plt.savefig('MF_v_SFE.png')
        print('made figure MF_v_SFE.png')
        
        plt.clf()
        plt.plot(SFE, Single_L_mean)
        plt.fill_between(SFE, np.array(Single_L_mean)-np.array(Single_L_std), np.array(Single_L_mean)+np.array(Single_L_std), alpha=0.2)
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.xlabel('SFE')
        plt.ylabel('Luminosity')
        plt.savefig('L_vs_SFE.png')
        
        plt.clf()
        plt.semilogy(SFE, Single_M_dot_mean)
        plt.fill_between(SFE, np.array(Single_M_dot_mean)-np.array(Single_M_dot_std), np.array(Single_M_dot_mean)+np.array(Single_M_dot_std), alpha=0.2)
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.xlabel('SFE')
        plt.ylabel('Accretion')
        plt.savefig('M_dot_vs_SFE.png')
        
sys.stdout.flush()
CW.Barrier()
#=====================================================
