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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-sim_G", "--simulation_G", type=str, default='')
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

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

sys.stdout.flush()
CW.Barrier()

#What do I need to save?
Times = []
SFE = []
MF = [[]]

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
            abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][n_stars])
            time = global_data['time'][time_it][n_stars][0]
            sfe = np.sum(mass)
            
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
            
            time_yt = yt.YTArray(time*scale_t, 's')
            Times.append(int(time_yt.in_units('yr').value))
            SFE.append(sfe)
            
            res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, Grho=Grho, max_iter=100)
            s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
            multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
            
            total_systems = len(s_true) + len(multi_inds)
            MF_value = len(multi_inds)/total_systems
            MF.append(MF_value)
            
            pickle_file_rank = pickle_file.split('.pkl')[0] + "_" +str(rank) + ".pkl"
            file = open(pickle_file_rank, 'wb')
            pickle.dump((Times, SFE, MF),file)
            file.close()
            print('updated pickle', pickle_file_rank, "for time_it", time_it, "of", end_time_ind+1)
            
sys.stdout.flush()
CW.Barrier()
print("FINISHED GOING THROUGH TIMES ON RANK", rank)
    
    
sys.stdout.flush()
CW.Barrier()
#=====================================================
