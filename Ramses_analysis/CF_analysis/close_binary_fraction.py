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
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-update", "--update_pickles", help="Do you want to remake the pickles?", type=str, default='True')
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


#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]
try:
    simulation_density_id_int = int(simulation_density_id)
except:
    simulation_density_id = args.global_data_pickle_file.split('_G')[-1].split('.pkl')[0]

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

#Define arrays of quantities to be saved.
SFE = []
Close_Fractions = []

sys.stdout.flush()
CW.Barrier()

#loading global data
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
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

#Get birth conditions
file_open = open(args.global_data_pickle_file, 'rb')
#Sink_bound_birth.append([born_bound, most_bound_sink_id, rel_sep])

sys.stdout.flush()
CW.Barrier()

if args.bound_check == 'True':
    bound_check = True
else:
    bound_check = False

sys.stdout.flush()
CW.Barrier()

#Find sink formation times
Sink_formation_times = []
time_it = 0
sink_it = 0
while time_it < np.shape(global_data['m'])[0]:
     while sink_it < np.shape(global_data['m'])[1]:
         if global_data['m'][time_it][sink_it]>0:
             Sink_formation_times.append(global_data['time'][time_it][sink_it])
             sink_it = sink_it + 1
         elif global_data['m'][time_it][sink_it] == 0:
             time_it = time_it + 1
             break
     if sink_it == np.shape(global_data['m'])[1]:
         break
Sink_formation_times = (Sink_formation_times*scale_t).in_units('yr')
            
sys.stdout.flush()
CW.Barrier()

start_time_it = args.start_time_index
plt.clf()
pickle_file = 'G'+simulation_density_id
if args.update_pickles == 'True':
    rit = -1
    prev_n_stars = 1
    for time_it in range(start_time_it, len(global_data['time'].T[0])):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            #preamble
            n_stars = np.where(global_data['m'][time_it]>0)[0]
            if len(n_stars) > 1:
                abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                mass = np.array(global_data['m'][time_it][n_stars])
                time = global_data['time'][time_it][n_stars][0]
                time_yr = yt.YTQuantity(scale_t*time, 's').in_units('yr').value
                
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
                
                M_tot_msun = np.sum(mass*units['mass_unit'].in_units('Msun'))
                SFE_val = M_tot_msun/units['mass_unit'].in_units('Msun')
                SFE.append(SFE_val)
                
                res = m.multipleAnalysis(S,cutoff=10000, bound_check=bound_check, nmax=6, cyclic=True, Grho=Grho, max_iter=args.max_iterations)
                
                multi_inds = np.where(res['n']>1)[0]
                close_seps = np.where(res['separation'][multi_inds] < 100)[0]
                if len(close_seps) == 0:
                    close_frac = 0
                else:
                    close_frac = len(close_seps)/len(multi_inds)
                    
                Close_Fractions.append(close_frac)
                
                print('calculated fration for time it', time_it, 'of', len(global_data['time'].T[0]))
                
                pickle_file_rank = pickle_file + '_' + str(rank) + '.pkl'
                file = open(pickle_file_rank, 'wb')
                pickle.dump((SFE, Close_Fractions),file)
                file.close()

sys.stdout.flush()
CW.Barrier()

#Compile pickles
if rank == 0:
    rank_pickles = glob.glob(pickle_file + '_*.pkl')
    SFE_all = []
    Close_Fractions_all = []
    for rank_pickle in rank_pickles:
        file = open(rank_pickle, 'rb')
        SFE, Close_Fractions = pickle.load(file)
        file.close()
        
        SFE_all = SFE_all + SFE
        Close_Fractions_all = Close_Fractions_all + Close_Fractions
        
    sort_inds = np.argsort(SFE_all)
    SFE = np.array(SFE_all)[sort_inds]
    Close_Fractions = np.array(Close_Fractions_all)[sort_inds]
        
    file = open(pickle_file+'.pkl', 'wb')
    pickle.dump((SFE, Close_Fractions),file)
    file.close()

sys.stdout.flush()
CW.Barrier()

plt.plot(SFE, Close_Fractions)
plt.ylim([0, 1])
plt.xlim([0, 0.05])
plt.xlabel('SFE')
plt.ylabel('Close binary fraction')
plt.savefig('close_binary_fraction_G'+simulation_density_id+'.png')
