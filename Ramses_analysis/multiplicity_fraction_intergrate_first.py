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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-proj", "--projected_separation", help="do you want to use projected separation instead of true separation?", default="False", type=str)
    parser.add_argument("-ax", "--axis", help="what axis do you want to project separation onto?", default='x', type=str)
    parser.add_argument("-preffix", "--figure_prefix", help="Do you want to give saved figures a preffix?", default="", type=str)
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/dev/shm/stars_red_512.pkl', type=str)
    parser.add_argument("-time_window", "--time_intergration_window", help="Over what time do you want to intergrate to calculate accretion rate?", default=1000.0, type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def losi(i, res):
    if (res['n'][i]==1):
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
        
#=====================================================================================================

#Baraffe models
L_bins = np.logspace(-1.55,1.55,11)
S_bins = np.logspace(1.25,4,12)
Tobin_luminosities_All_objects = np.array([0.04,0.05,0.09,0.1,0.1,0.16,0.16,0.16,0.17,0.23,0.24,0.25,0.3,0.3,0.3,0.32,0.36,0.38,0.39,0.4,0.4,0.43,0.5,0.5,0.54,0.54,0.54,0.54,0.6,0.6,0.63,0.68,0.69,0.7,0.7,0.7,0.8,0.87,0.9,1,1.1,1.2,1.2,1.3,1.3,1.4,1.4,1.4,1.5,1.5,1.5,1.6,1.7,1.8,1.8,1.8,1.8,1.9,16.8,19,2.5,2.6,2.8,23.2,3.2,3.2,3.6,3.7,32.5,4,4.2,4.7,5.3,6.9,7,8.3,8.4,9.1,9.2,0.04,0.04,0.04,0.05,0.05,0.05,0.05,0.07,0.07,0.14,0.15,0.22,0.28,0.47,1.1,1.3])
Tobin_Luminosities_multiples = np.array([1.02,1.4,0.04,3.2,1.3,1.4,1.5,1.5,7,11,4.1,1.1,4.2,2.8,3.9,0.9,9.7,3.6,9.08,19,0.3,2.1,8.3,8.3,17.5,9.7,9.1,5.3,24.3,1.9,
1.04,1.5,32.5,33.5,34,0.87,1.1,1.3,1.8,0.79,0.9,4.4])
Tobin_hist, bins = np.histogram(Tobin_Luminosities_multiples, bins=L_bins)
        
#=====================================================================================================

args = parse_inputs()

units = {"length_unit":yt.YTQuantity(4.0,"pc"), "mass_unit":yt.YTQuantity(2998,"Msun"), "velocity_unit":yt.YTQuantity(0.18, "km/s"), "time_unit":yt.YTQuantity(685706129102738.9, "s"), "density_unit":yt.YTQuantity(46.84375, "Msun/pc**3")}

scale_l = 1.23427103e19 # 4 pc
scale_v = 1.8e4         # 0.18 km/s == sound speed
scale_t = 6.85706128e14 # 4 pc / 0.18 km/s
scale_d = 3.171441e-21  # 2998 Msun / (4 pc)^3

window = yt.YTQuantity(args.time_intergration_window, 'yr')

'''
#low resolution data
datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data_256'
nout = 133
'''

#high resolution data
#datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data'
datadir = sys.argv[1]
savedir = sys.argv[2]
file_no_range = [298,407]

#Load global pickle data
pickle_file = args.global_data_pickle_file
file_open = open(pickle_file, 'rb')
global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
print('Loaded global pickle data')

#Get total number of sinks:
S = pr.Sink()
S.nout = file_no_range[-1]
S.datadir = datadir
S.rsink()
S._jet_factor = 1.
S._scale_l = scale_l
S._scale_v = scale_v
S._scale_t = scale_t
S._scale_d = scale_d
i = S.read_info_file(file_no_range[-1],datadir=datadir)
S._time = i['time']
nmax=len(S.mass)

#Define quantities to keep track off
Separations = []
Times = []
Luminosities = []
N_sys_total = []
CF_Array_Full = []
All_unique_systems = {}
All_unique_systems_L = {}
All_unique_systems_M = {}

for nout in range(file_no_range[0], file_no_range[1]+1):
    L_tots = []
    usable_inds = []

    # load sink data from snapshot nout in to S
    S = pr.Sink()
    S.nout = nout
    S.datadir = datadir
    #Get data from previous time step and next time step:
    S.rsink()
    S._jet_factor = 1.
    S._scale_l = scale_l
    S._scale_v = scale_v
    S._scale_t = scale_t
    S._scale_d = scale_d
    i = S.read_info_file(nout,datadir=datadir)
    S._time = i['time']
    res = m.multipleAnalysis(S,nmax=nmax,cutoff=1e4)
    time = S._time * units['time_unit'].in_units('yr')
    top_systems = np.where(res['topSystem'])[0]
    
    s = np.where(res['n'][top_systems]==1)[0]
    b = np.where(res['n']==2)[0]
    t = np.where(res['n']==3)[0]
    q = np.where(res['n']==4)[0]
    q5 = np.where(res['n']==5)[0]
    s6 = np.where(res['n']==6)[0]
    h = np.where(res['n']==7)[0]
    o = np.where(res['n']==8)[0]
    higher_order = np.where(res['n']>8)[0]
    multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist() + h.tolist() + o.tolist() + higher_order.tolist())

    #CALCULATE MAGNITUDE
    radius = yt.YTQuantity(2.0, 'rsun')
    temperature = yt.YTQuantity(3000, 'K')
    L_phot = yt.units.stefan_boltzmann_constant * 4*np.pi*radius.in_units('cm')**2 * temperature**4
    
    #Find Luminosity of system
    min_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time+window))
    max_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time-window))
    global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time))
    
    #Find Single stars that are withink the limonsity limits. Single stars must be Top Systems.
    
    for s_ind in s:
        dM = (global_data['m'][max_time_ind,s_ind] - global_data['m'][min_time_ind,s_ind])*units['mass_unit'].in_units('msun')
        dt = (global_data['time'][max_time_ind,s_ind] - global_data['time'][min_time_ind,s_ind])*units['time_unit'].in_units('yr')
        M_dot = dM/dt
        M = yt.YTArray(global_data['m'][global_ind,s_ind]*units['mass_unit'].in_units('msun'), 'Msun')
        f_acc = 0.5
        L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
        L_tot = L_acc.in_units('Lsun')
        mean_L = np.sum(L_tot)
        
        if str([s_ind]) not in All_unique_systems.keys():
            All_unique_systems.update({str([s_ind]): [res['separation'][s_ind]]})
            All_unique_systems_M.update({str([s_ind]): [M]})
            All_unique_systems_L.update({str([s_ind]): mean_L})
        else:
            All_unique_systems[str([s_ind])].append(res['separation'][s_ind])
            All_unique_systems_M[str([s_ind])].append(M)
            All_unique_systems_L[str([s_ind])] = np.append(All_unique_systems_L[str([s_ind])],mean_L)
        
    
    #Then add the data for the multiple systems. They don't necessarily need to be Top Systems
    
    for multi_ind in multi_inds:
        sys_comps = losi(multi_ind, res)
        sys_comps = sorted(flatten(sys_comps))
        dM = (global_data['m'][max_time_ind,sys_comps] - global_data['m'][min_time_ind,sys_comps])*units['mass_unit'].in_units('msun')
        dt = (global_data['time'][max_time_ind,sys_comps] - global_data['time'][min_time_ind,sys_comps])*units['time_unit'].in_units('yr')
        M_dot = dM/dt
        M = yt.YTArray(global_data['m'][global_ind,sys_comps]*units['mass_unit'].in_units('msun'), 'Msun')
        f_acc = 0.5
        L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
        L_tot = L_acc.in_units('Lsun')
        mean_L = np.sum(L_tot)
        
        if str(sys_comps) not in All_unique_systems.keys():
            All_unique_systems.update({str(sys_comps): [res['separation'][multi_ind]]})
            All_unique_systems_M.update({str(sys_comps): [M]})
            All_unique_systems_L.update({str(sys_comps): mean_L})
        else:
            All_unique_systems[str(sys_comps)].append(res['separation'][multi_ind])
            All_unique_systems_M[str(sys_comps)].append(M)
            All_unique_systems_L[str(sys_comps)] = np.append(All_unique_systems_L[str(sys_comps)],mean_L)
            
    print("added data for nout=", nout)
        
   
