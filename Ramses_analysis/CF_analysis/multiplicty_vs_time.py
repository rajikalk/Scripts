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
Single_L_min = []
Single_L_max = []
Single_M_dot_min = []
Single_M_dot_max = []

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
    
SFE = np.sum(global_data['m'], axis=1)
SFE_vals = [0.01, 0.02, 0.03, 0.04, 0.05]
SFE_window = 0.001
time_its = []
for SFE_val in SFE_vals:
    start_SFE = SFE_val - SFE_window
    end_SFE = SFE_val + SFE_window
    start_time_it = np.argmin(abs(SFE-start_SFE))
    end_time_it = np.argmin(abs(SFE-end_SFE))
    time_its = np.arange(start_time_it, end_time_it)
    for time_it in time_its:
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
         
        import pdb
        pdb.set_trace()
