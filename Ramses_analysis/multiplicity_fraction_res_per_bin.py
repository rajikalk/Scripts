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
f_acc = 0.5

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-proj", "--projected_separation", help="do you want to use projected separation instead of true separation?", default="False", type=str)
    parser.add_argument("-ax", "--axis", help="what axis do you want to project separation onto?", default='x', type=str)
    parser.add_argument("-preffix", "--figure_prefix", help="Do you want to give saved figures a preffix?", default="", type=str)
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/stars_red_512.pkl', type=str)
    parser.add_argument("-time_window", "--time_intergration_window", help="Over what time do you want to intergrate to calculate accretion rate?", default=1000.0, type=float)
    parser.add_argument("-verbose", "--verbose_printing", help="Would you like to print debug lines?", type=str, default='False')
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=35.6)
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-lifetime", "--lifetime_threshold", help="What life time threshold do you want to consider when making Luminosity histogram", type=float, default=10000)
    parser.add_argument("-proj_vec", "--projection_vector", help="What projection vector do you want to use?", type=str, default='')
    parser.add_argument("-plot_only", "--make_plots_only", help="Do you just want to make plots? Not calculate the CF", type=str, default='False')
    parser.add_argument("-t_spread", "--time_spread", help="how much time around the central time do you want to intergrate over?", type=float, default=0.01)
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

def accretion(global_data, sink_inds, max_time_ind, min_time_ind):
    """
    Calculates the accretion of the given indeexes
    """
    dM = (global_data['m'][max_time_ind,sink_inds] - global_data['m'][min_time_ind,sink_inds])*units['mass_unit'].in_units('msun')
    dt = (global_data['time'][max_time_ind,sink_inds] - global_data['time'][min_time_ind,sink_inds])*units['time_unit'].in_units('yr')
    M_dot = dM/dt
    return M_dot
    

def luminosity(global_data, sink_inds, max_time_ind, min_time_ind):
    """
    Calculates the luminosity of the given indexes
    """
    global f_acc
    M_dot = accretion(global_data, sink_inds, max_time_ind, min_time_ind)
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
cpus_per_time = len(S_bins) - 1
#S_bins = np.logspace(1,4,13)
#S_bins = np.logspace(1.25,4,12)
Tobin_luminosities_All_objects = np.array([0.04,0.05,0.09,0.1,0.1,0.16,0.16,0.16,0.17,0.23,0.24,0.25,0.3,0.3,0.3,0.32,0.36,0.38,0.39,0.4,0.4,0.43,0.5,0.5,0.54,0.54,0.54,0.54,0.6,0.6,0.63,0.68,0.69,0.7,0.7,0.7,0.8,0.87,0.9,1,1.1,1.2,1.2,1.3,1.3,1.4,1.4,1.4,1.5,1.5,1.5,1.6,1.7,1.8,1.8,1.8,1.8,1.9,16.8,19,2.5,2.6,2.8,23.2,3.2,3.2,3.6,3.7,32.5,4,4.2,4.7,5.3,6.9,7,8.3,8.4,9.1,9.2,0.04,0.04,0.04,0.05,0.05,0.05,0.05,0.07,0.07,0.14,0.15,0.22,0.28,0.47,1.1,1.3])
Tobin_Luminosities_multiples = np.array([0.9, 1.3, 4.2, 0.87, 1.5, 1.3, 3.6, 3.2, 9.1, 0.04, 9.08, 1.5, 4.4, 1.1, 1.19, 18.9, 10.8, 0.79, 11.1, 24.3, 0.9, 2.44, 35.54, 2.1])

Tobin_objects = {
"Per-emb-2": [18.4],
"Per-emb-5": [22.3],
"Per-emb-17": [63.9],
"Per-emb-48": [79.5],
"Per-emb-40": [90],
"EDJ2009-269": [120.6],
"Per-emb-22": [172.8],
"EDJ2009-183": [235.8],
"L1448IRS1": [327.4],
"Per-emb-35": [438.8],
"EDJ2009-156": [714.6],
"Per-emb-26+Per-emb-42": [1864],
"Per-emb-11": [678.8, 2177.8],
"Per-emb-8+Per-emb-55": [142.1, 2198.2],
"Per-emb-16+Per-emb-28": [3694.5],
"B1-bN+B1-bS+Per-emb-41": [3210.1, 4000.8],
"Per-emb-33+L1448NW+L1448IRS3A": [57.7, 60.7, 182.8, 1683, 4945.6],
"Per-emb-21+Per-emb-18+Per-emb-49": [19.6, 71.9, 3048, 9367.1],
"Per-emb-58+Per-emb-65": [6641.9],
"Per-emb-12+Per-emb-13+IRAS4B'": [420.8, 2450.4, 6840],
"Per-emb-36+Per-emb-27": [71.6, 142.6, 7226.6],
"Per-emb-6+Per-emb-10": [7347.9],
"Per-emb-37+EDJ2009+233+EDJ2009+235": [2427.8, 7752],
"Per-emb-44+SVS13A2+SVS13B+SVS13C": [69, 1222.2, 3434.4, 7941.5],
"Per-emb-32+EDJ2009+366": [1395.3, 8419.2]}

Tobin_hist, bins = np.histogram(Tobin_Luminosities_multiples, bins=L_bins)

MF_Tobin = []
Tobin_intersection = list(set(Tobin_luminosities_All_objects).intersection(Tobin_Luminosities_multiples))
Tobin_n_singles = 60
N_components = np.ones(Tobin_n_singles).tolist()
CF_per_bin_Tobin = []
CF_errs = []

sys.stdout.flush()
CW.Barrier()

for bin_it in range(1,len(S_bins)):
    N_components = []#np.ones(Tobin_n_singles).tolist()
    for key in Tobin_objects.keys():
        bin_inds = np.where((np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]))[0]
        if len(np.where(Tobin_objects[key]<S_bins[bin_it])[0]) == 0:
            #If all separations are greater than bin upper bound, every component is taken to be a single star
            N_components = N_components + np.ones(len(Tobin_objects[key])+1).tolist()
        elif len(np.where((np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]))[0]) == len(Tobin_objects[key]):
            #If all components are in the separation bin, then N = number of components
            N_components.append(len(Tobin_objects[key])+1)
        elif len(np.where(Tobin_objects[key]<S_bins[bin_it-1])[0]) == 0:
            #If all separations are smaller than the bin upper bound then the whole system is taken to be a single star.
            N_components.append(1)
        elif True not in (np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]):
            #If some separations are either above or below the bin bounds, then separations above the upper bin bound are single stars, plus one for the separations below the lower bound.
            n_single_comps = len(np.where(Tobin_objects[key]>S_bins[bin_it])[0]) + 1
            N_components = N_components + np.ones(n_single_comps).tolist()
        elif bin_inds[0] != 0 and bin_inds[-1] != len(Tobin_objects[key])-1:
            #if some of the middle separations are in the bin, we calculate what would be seen as singles and what is the oerder of the observed multiple would be
            bin_inds = np.where((np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]))[0]
            n_single_comps = len(Tobin_objects[key][bin_inds[-1]+1:])
            N_components = N_components + np.ones(n_single_comps).tolist()
            N_components.append(len(bin_inds)+1)
        elif ((np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]))[0] == False and ((np.array(Tobin_objects[key])<S_bins[bin_it])&(np.array(Tobin_objects[key])>S_bins[bin_it-1]))[-1] == True:
            N_components.append(len(bin_inds)+1)
        else:
            import pdb
            pdb.set_trace()
    s_no = len(np.where(np.array(N_components) == 1)[0])
    b_no = len(np.where(np.array(N_components) == 2)[0])
    t_no = len(np.where(np.array(N_components) == 3)[0])
    q_no = len(np.where(np.array(N_components) == 4)[0])
    cf = (b_no+t_no*2+q_no*3)/(s_no+b_no+t_no+q_no)
    CF_per_bin_Tobin.append(cf)
    
    N_comp = np.sum(np.array(N_components) - np.ones(np.shape(N_components)))
    N_sys = len(N_components)
    CF_err = (1/N_sys)*(N_comp*(1-(N_comp/N_sys)))**0.5
    CF_errs.append(CF_err)

sys.stdout.flush()
CW.Barrier()

CF_per_bin_Tobin = np.array(CF_per_bin_Tobin)
CF_errs = np.array(CF_errs)

#Create histogram of Tobin's CF Distribution
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
#plt.xlim([1,4])
#plt.ylim([0, 0.2])
plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
plt.ylim(bottom=0)

#raghaven dist
sep_mean_rag = 1.7
sep_std_rag = 1.52
x = np.linspace(0.75, 4, 1000)
CF_rag = norm.pdf(x,sep_mean_rag,sep_std_rag)*0.15240300488036065
#plt.plot(x, CF_rag)
CF_rag_des = []
for bin_it in range(1,len(S_bins)):
    bin_inds = np.where(((10**x)<S_bins[bin_it])&((10**x)>S_bins[bin_it-1]))[0]
    CF_bin_rag = np.mean(CF_rag[bin_inds])
    CF_rag_des.append(CF_bin_rag)

def Gaussian(x,scale,mean,sigma):
    return scale * norm.pdf(x, mean, sigma)

diff_hist = np.array(CF_per_bin_Tobin) - np.array(CF_rag_des)
diff_hist = np.clip(diff_hist, 0, np.max(diff_hist))
#plt.bar(bin_centers, diff_hist, yerr=CF_errs, width=0.25, edgecolor='black', alpha=0.5)
popt, pcov = curve_fit(Gaussian, bin_centers[-5:], diff_hist[-5:], [0.3, 4.0, 0.3], sigma=CF_errs[-5:])

fit = Gaussian(x, *popt)

#Tobin Gaussian fits
first_peak = np.log10(90)
popt1, pcov1 = curve_fit(Gaussian, bin_centers[:9], CF_per_bin_Tobin[:9], [np.max(CF_per_bin_Tobin[:9]), bin_centers[:9][np.argmax(CF_per_bin_Tobin[:9])], 0.5], sigma=CF_errs[:9])
fit1 = Gaussian(x, *popt1)
plt.plot(x, fit1, ls='--')

start_ind = 8
popt2, pcov2 = curve_fit(Gaussian, bin_centers[start_ind:], CF_per_bin_Tobin[start_ind:], [np.max(CF_per_bin_Tobin[start_ind:]), bin_centers[start_ind:][np.argmax(CF_per_bin_Tobin[start_ind:])], 0.5], sigma=CF_errs[start_ind:])
fit2 = Gaussian(x, *popt2)
plt.plot(x, fit2, ls='--')

sum_fit = fit1 + fit2
plt.plot(x, sum_fit)

def line(x, m, b):
    return m*x + b


non_zero_inds = np.where(CF_per_bin_Tobin>0)[0]
popt3, pcov3 = curve_fit(line, bin_centers[non_zero_inds][:7], CF_per_bin_Tobin[non_zero_inds][:7], [0, np.mean(CF_per_bin_Tobin[non_zero_inds][:7])], sigma=CF_errs[non_zero_inds][:7])
popt3[0] = 0
popt3[1] = np.mean(CF_per_bin_Tobin[non_zero_inds][:7])
fit3 = line(x, *popt3)
plt.plot(x, fit3, ls='--')
popt4, pcov4 = curve_fit(line, bin_centers[9:], CF_per_bin_Tobin[9:]-popt3[1], [ 0.39, -1.15], sigma=CF_errs[9:])
#power_lw = power_law(x, 0.00004, 10, 0.025)
fit4 = line(x, *popt4)
plt.plot(x, fit4, ls='--')

fit4 = np.clip(fit4, 0, np.max(fit4))
sum_line = fit3 + fit4
plt.plot(x, sum_line)
#sep_mean = 4.0
#sep_std = 0.3
#plt.plot(x, norm.pdf(x,sep_mean,sep_std)*0.3)

plt.savefig(savedir + "CF_Tobin_data.png")

sys.stdout.flush()
CW.Barrier()
        
#file = open("Tobin_CF.pkl", 'wb')
#pickle.dump((S_bins, CF_per_bin_Tobin),file)
#print('updated pickle', "Tobin_CF.pkl")
#file.close()

#=====================================================================================================

if len(args.projection_vector) > 0:
    proj_vector = np.array(eval(args.projection_vector))
    proj_length = np.sqrt(np.sum(proj_vector**2))
    proj_unit = proj_vector/proj_length
else:
    proj_unit = []

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.key():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

sys.stdout.flush()
CW.Barrier()

if args.pickled_file != None:
    pickle_file = savedir + args.pickled_file
    if os.path.isfile(pickle_file):
        file = open(pickle_file, 'rb')
        Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
        file.close()
    else:
        #Define quantities to keep track off
        Separations = []
        Times = []
        Luminosities = []
        N_sys_total = []
        CF_Array_Full = []
        All_unique_systems = {}
        All_unique_systems_L = {}
        All_unique_systems_M = {}
        All_unique_systems_T = {}
else:
    #Define quantities to keep track off
    Separations = []
    Times = []
    Luminosities = []
    N_sys_total = []
    CF_Array_Full = []
    All_unique_systems = {}
    All_unique_systems_L = {}
    All_unique_systems_M = {}
    All_unique_systems_T = {}
    
sys.stdout.flush()
CW.Barrier()

if args.bound_check == 'True':
    bound_check = True
else:
    bound_check = False

#Calculate variables
window = yt.YTQuantity(args.time_intergration_window, 'yr')
luminosity_lower_limit = 0.01
luminosity_upper_limit = args.upper_L_limit
accretion_limit = args.accretion_limit

'''
#low resolution data
datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data_256'
nout = 133
'''

#high resolution data
#datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data'

#Load global pickle data
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

SFE_ind = np.argmin(np.abs(0.049-(np.sum(global_data['m'], axis=1)*units['mass_unit'].value)/units['mass_unit'].value))
SFE_t_ff = global_data['time'].T[0][SFE_ind]
dt = args.time_spread
time_bounds = [SFE_t_ff-dt,SFE_t_ff+dt]

start_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[0]))
end_time_ind = np.argmin(abs(global_data['time'].T[0]-time_bounds[1]))

print("tstart, tend", time_bounds[0]*units['time_unit'].in_units('kyr'), time_bounds[1]*units['time_unit'].in_units('kyr'))
print("SFE_start, SFE_end", (np.sum(global_data['m'][start_time_ind])*units['mass_unit'].value)/units['mass_unit'].value, (np.sum(global_data['m'][end_time_ind])*units['mass_unit'].value)/units['mass_unit'].value)
print("nstars_start, nstars_end", len(np.where(global_data['m'][start_time_ind]>0)[0]), len(np.where(global_data['m'][end_time_ind]>0)[0]))
print("mstars_start, mstars_end", np.sum(global_data['m'][start_time_ind])*units['mass_unit'].value, np.sum(global_data['m'][end_time_ind])*units['mass_unit'].value)
print("number of records", end_time_ind-start_time_ind)

radius = yt.YTQuantity(2.0, 'rsun')
temperature = yt.YTQuantity(3000, 'K')

if len(Times) == (end_time_ind-start_time_ind):
    update = False
else:
    update = True
    start_time_ind = start_time_ind + len(Times)
    
sys.stdout.flush()
CW.Barrier()

if update == True and args.make_plots_only == 'False':
    for time_it in range(start_time_ind, end_time_ind+1):
        L_tots = []
        current_separations = []

        # Make S with the data at time time_it
        n_stars = np.where(global_data['m'][time_it]>0)[0]
        abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
        absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
        mass = np.array(global_data['m'][time_it][n_stars])
        time = global_data['time'][time_it][n_stars][0]

        S = pr.Sink()
        S._jet_factor = 1.
        S._scale_l = scale_l
        S._scale_v = scale_v
        S._scale_t = scale_t
        S._scale_d = scale_d
        S._time = yt.YTArray(time, '')
        S._abspos = yt.YTArray(abspos, '')
        S._absvel = yt.YTArray(absvel, '')
        S._mass = yt.YTArray(mass, '')

        time_yt = yt.YTArray(time*scale_t, 's')
        Times.append(time_yt.in_units('yr').value)
        
        #CALCULATE Photospheric Luminosity
        L_phot = yt.units.stefan_boltzmann_constant * 4*np.pi*radius.in_units('cm')**2 * temperature**4
        
        #Find indexes of time window
        min_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')+window))
        max_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')-window))
        global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time_yt.in_units('yr')))
        
        #Save cf and system breakdown (S:B:T:Q:..) for each bin
        CF_per_bin = []
        n_systems = []
        
        sys.stdout.flush()
        CW.Barrier()
        
        #Iterate over each bin:
        if size == 1:
            rit = 0
        else:
            rit = -1
        for bin_it in range(1,len(S_bins)):
            if size != 1:
                rit = rit + 1
            if rank == rit:
                sep_update = {}
                L_update = {}
                T_update = {}
                #Calculate res for each bin upper bound
                if args.projected_separation == "False":
                    #res = m.multipleAnalysis(S,nmax=6,cutoff=S_bins[bin_it], bound_check=bound_check, Grho=Grho)
                    res = m.multipleAnalysis(S,cutoff=S_bins[bin_it], bound_check=bound_check, nmax=6, cyclic=False, Grho=Grho)
                else:
                    #res = m.multipleAnalysis(S,nmax=6,cutoff=S_bins[bin_it], bound_check=bound_check, projection=True, axis=args.axis, projection_vector=proj_unit, Grho=Grho)
                    res = m.multipleAnalysis(S,cutoff=S_bins[bin_it], bound_check=bound_check, nmax=6, cyclic=False, projection=True, axis=args.axis, projection_vector=proj_unit, Grho=Grho)
                if size == 1:
                    if np.max(res['n']) > 6:
                        import pdb
                        pdb.set_trace()
                top_inds = np.where(res['topSystem'])[0]

                sink_inds = np.where((res['n']==1))[0]
                L_tot = luminosity(global_data, sink_inds, max_time_ind, min_time_ind)
                M_dot = accretion(global_data, sink_inds, max_time_ind, min_time_ind)
                s_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<luminosity_upper_limit))[0]
                visible_stars = sink_inds[s_inds]
                
                if args.verbose_printing != 'False':
                    print("=========================================================")
                    print("TRUE NUMBER OF VISIBLE STARS IS", str(len(visible_stars)))

                #Find all singles and top systems with separations below the bin lower bound
                s_true = np.where((res['n']==1) & (res['topSystem']==True))[0]
                s_fake = np.where((res['separation']<S_bins[bin_it-1])&(res['topSystem']==True)&(res['n']!=1))[0]
                
                #Filter out invisible systems
                #Find Single stars that are within the limonsity limits
                L_tot = luminosity(global_data, s_true, max_time_ind, min_time_ind)
                M_dot = accretion(global_data, s_true, max_time_ind, min_time_ind)
                s_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<luminosity_upper_limit))[0]#&(M_dot>accretion_limit))[0] #&(L_tot<35.0)
                visible_stars = s_true[s_inds]
                if args.verbose_printing != 'False':
                    print("AND", len(visible_stars), "ARE SINGLE STARS")
                invisible_stars = list(set(s_true).symmetric_difference(visible_stars))
                if args.verbose_printing != 'False':
                    print("UP TO " + str(int(S_bins[bin_it])) +"AU:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds])))
                res['n'][invisible_stars] = 0
                if args.verbose_printing != 'False':
                    print("AFTER REMOVING " + str(len(invisible_stars)) + " INVISIBLE SINGLE STARS:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds])))
                
                removed_sys = 0
                redefined_n = 0
                stars_in_collapsed_systems = 0
                for fake in s_fake:
                    sys_comps = losi(fake, res)
                    sys_comps = sorted(flatten(sys_comps))
                    
                    #Find luminosities of individual components
                    L_tot = luminosity(global_data, sys_comps, max_time_ind, min_time_ind)
                    M_dot = accretion(global_data, sys_comps, max_time_ind, min_time_ind)
                    visible_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<luminosity_upper_limit))[0] #&((dM/dt)>accretion_limit))[0]
                    visible_components = np.array(sys_comps)[visible_inds]
                    invisible_stars = list(set(sys_comps).symmetric_difference(visible_components))
                    
                    res['n'][visible_components] = 1
                    res['n'][invisible_stars] = 0
                    
                    if len(visible_components) > 0: #&(M_dot>accretion_limit):
                        stars_in_collapsed_systems = stars_in_collapsed_systems + len(visible_components)
                        res['n'][fake] = 1
                        redefined_n = redefined_n + 1
                    else:
                        res['n'][fake] = 0
                        removed_sys = removed_sys + 1
                
                missing_stars = stars_in_collapsed_systems - redefined_n
                if args.verbose_printing != 'False':
                    print("AFTER REMOVING " + str(removed_sys) + " SYSTEMS, REFINING " + str(redefined_n) + " SYSTEMS (WITH " + str(stars_in_collapsed_systems) + " STARS):" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds])))
                del removed_sys
                del redefined_n
                del stars_in_collapsed_systems
                del s_fake
                del s_true
                
                #Determine which systems could still be multiples
                b = np.where((res['n']==2) & (res['topSystem']==True))[0]
                t = np.where((res['n']==3) & (res['topSystem']==True))[0]
                q = np.where((res['n']==4) & (res['topSystem']==True))[0]
                q5 = np.where((res['n']==5) & (res['topSystem']==True))[0]
                s6 = np.where((res['n']==6) & (res['topSystem']==True))[0]
                s7 = np.where((res['n']==7) & (res['topSystem']==True))[0]
                multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist() + s7.tolist())

                del b
                del t
                del q
                del q5
                del s6
                del s7
            
                #Now go over multiple systems
                for multi_ind in multi_inds:
                    sys_comps = losi(multi_ind, res)
                    sys_comps = sorted(flatten(sys_comps))
                    
                    L_tot = luminosity(global_data, sys_comps, max_time_ind, min_time_ind)
                    M_dot = accretion(global_data, sys_comps, max_time_ind, min_time_ind)
                    
                    detectable_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit)&(L_tot<luminosity_upper_limit))[0]
                    detectable_components = np.array(sys_comps)[detectable_inds]
                                
                    #Recalculate luminsoity and accretion for detectable inds
                    M_dot = accretion(global_data, detectable_components, max_time_ind, min_time_ind)
                    L_tot = luminosity(global_data, detectable_components, max_time_ind, min_time_ind)
                    mean_M_dot = np.sum(M_dot)
                    mean_L = np.sum(L_tot)
                    
                    #if (mean_L+std_L) >0.04 and (mean_L-std_L) < 32.0 and mean_M_dot+std_M_dot>1.e-7:
                    if mean_L >luminosity_lower_limit and mean_M_dot>accretion_limit and mean_L<luminosity_upper_limit:
                        #Check whether there is undetected companions
                        n_detectable = len(detectable_components)
                        res['n'][multi_ind] = n_detectable
                        if n_detectable == len(sys_comps):
                            #Check whether hierarchial. If it is it should be pretty easy for code:
                            system_structure = losi(multi_ind, res)
                        else:
                            system_structure = losi(multi_ind, res)
                            invisible_components = list(set(sys_comps).symmetric_difference(detectable_components))
                            res['n'][invisible_components] = 0
                            if n_detectable == 0:
                                print("none of the components are detected")
                            elif n_detectable == 1:
                                sys_comps = detectable_components
                            elif n_detectable == 2:
                                system_structure = detectable_components
                                sys_comps = sorted(flatten(system_structure))
                            else:
                                sys_string = str(system_structure)
                                for inv_comp in invisible_components:
                                    inv_string = '['+str(inv_comp)+','
                                    if len(sys_string.split('['+str(inv_comp)+',')) > 1:
                                        split_string = sys_string.split(inv_string)
                                        sys_string = '[,'.join(split_string)
                                    else:
                                        inv_string = ' '+str(inv_comp)+']'
                                        split_string = sys_string.split(inv_string)
                                        sys_string = ' ]'.join(split_string)
                                
                                reduced = False
                                while reduced == False:
                                    if ', []' in sys_string:
                                        sys_string = ''.join(sys_string.split(', []'))
                                    if '[], ' in sys_string:
                                        sys_string = ''.join(sys_string.split('[], '))
                                    reduced_brackets = False
                                    while reduced_brackets == False:
                                        open_braket_inds = []
                                        delete_bool = []
                                        for char_it in range(len(sys_string)):
                                            if sys_string[char_it] == '[':
                                                open_braket_inds.append(char_it)
                                                if sys_string[char_it:char_it+3] == '[, ':
                                                    delete_bool.append(True)
                                                else:
                                                    delete_bool.append(False)
                                            elif sys_string[char_it] == ']':
                                                open_ind = open_braket_inds.pop()
                                                del_bool = delete_bool.pop()
                                                if sys_string[char_it-2:char_it+1] == ', ]':
                                                    str_1 = sys_string[:open_ind]
                                                    str_2 = sys_string[open_ind+1:char_it-2]
                                                    str_3 = sys_string[char_it+1:]
                                                    sys_string = str_1 + str_2 + str_3
                                                    break
                                                elif del_bool == True:
                                                    str_1 = sys_string[:open_ind]
                                                    str_2 = sys_string[open_ind+3:char_it]
                                                    str_3 = sys_string[char_it+1:]
                                                    sys_string = str_1 + str_2 + str_3
                                                    break
                                        if False not in delete_bool:
                                            reduced_brackets = True
                                    try:
                                        system_structure = eval(sys_string)
                                        reduced = True
                                    except:
                                        print("sys_string", sys_string, "is still not reduced")
                                
                                sys_comps = sorted(flatten(system_structure))
                                        
                        #Now calculate separations according to Tobin
                        if n_detectable == 1 or n_detectable == 0:
                            sep_list = []
                        
                        if n_detectable == 1:
                            L_list = [L_tot[0]]
                        elif n_detectable == 0:
                            L_list = []
                        elif is_hierarchical(system_structure):
                            central_ind = sys_comps[np.argmax(L_tot)]
                            positions = res['abspos'][sys_comps]
                            central_pos = res['abspos'][central_ind]
                            separations = np.sqrt(np.sum((positions - central_pos)**2, axis=1))
                            non_zero_inds = np.where(separations>0)[0]
                            sep_list = sorted(separations[non_zero_inds].tolist())
                            proximity_inds = np.argsort(separations)
                            p_ind = 2
                            L_list = [L_tot[proximity_inds[0]] + L_tot[proximity_inds[1]]]
                            while p_ind < len(proximity_inds):
                                L_list.append(L_list[-1] + L_tot[proximity_inds[p_ind]])
                                p_ind = p_ind + 1

                            invisible_sub_systems = np.where(sep_list < S_bins[bin_it-1])[0]
                            if len(invisible_sub_systems) > 0:
                                #print("Systems is hierarchical and there is a separation that is below the bin lower limit")
                                
                                separations_to_remove = np.array(sep_list)[invisible_sub_systems]
                                for sep_rm in separations_to_remove:
                                    sep_list.remove(sep_rm)
                                    
                                L_list = L_list[-1*len(sep_list):]
                        else:
                            sep_list = []
                            L_list = []
                            if len(flatten(system_structure)) == 4 and np.shape(system_structure) == (2,2):
                                primary_positions = []
                                for sys_s in system_structure:
                                    binary_positions = res['abspos'][sys_s]
                                    binary_sep = np.sqrt(np.sum((binary_positions[0] -  binary_positions[1])**2))
                                    sep_list.append(binary_sep)
                                    sys_comps_ind = []
                                    for comp in sys_s:
                                        sys_comps_ind.append(sys_comps.index(comp))
                                    primary_ind = sys_s[np.argmax(np.array(L_tot)[sys_comps_ind])]
                                    primary_positions.append(res['abspos'][primary_ind])
                                    L_list.append(np.sum(L_tot[sys_comps_ind]))
                                quad_sep = np.sqrt(np.sum((primary_positions[0] -  primary_positions[1])**2))
                                sep_list.append(quad_sep)
                                L_list.append(np.sum(L_tot))
                                
                                invisible_sub_systems = np.where(sep_list < S_bins[bin_it-1])[0]
                                if len(invisible_sub_systems) > 0:
                                    separations_to_remove = np.array(sep_list)[invisible_sub_systems]
                                    for sep_rm in separations_to_remove:
                                        sep_list.remove(sep_rm)
                                    
                                    luminosity_to_remove = np.array(L_list)[invisible_sub_systems]
                                    for L_rm in luminosity_to_remove:
                                        L_list.remove(L_rm)
                                    
                            else:
                                sep_list = []
                                L_list = []
                                raw_sys_string = str(system_structure)
                                sys_string = str(system_structure)
                                open_braket_inds = []
                                reduced = False
                                while reduced == False:
                                    for char_it in range(len(raw_sys_string)):
                                        if raw_sys_string[char_it] == '[':
                                            open_braket_inds.append(char_it)
                                        elif raw_sys_string[char_it] == ']':
                                            open_ind = open_braket_inds.pop()
                                            #found a sub-system, now reduce
                                            try:
                                                sub_system = eval(raw_sys_string[open_ind:char_it+1])
                                            except:
                                                sub_system = eval(raw_sys_string[open_ind:char_it])
                                            if len(sub_system) == 2:
                                                binary_positions = res['abspos'][sub_system]
                                                binary_sep = np.sqrt(np.sum((binary_positions[0] - binary_positions[1])**2))
                                                sep_list.append(binary_sep)
                                                sys_comps_ind = []
                                                for comp in sub_system:
                                                    sys_comps_ind.append(sys_comps.index(comp))
                                                primary_ind = sub_system[np.argmax(np.array(L_tot)[sys_comps_ind])]
                                                L_list.append(np.sum(L_tot[sys_comps_ind]))
                                                sys_string = str(primary_ind).join(sys_string.split(str(sub_system)))
                                                try:
                                                    reduced_system = eval(sys_string)
                                                    if is_hierarchical(reduced_system):
                                                        reduced = True
                                                        break
                                                    else:
                                                        print("sys_string", sys_string, "is still non-hierachical")
                                                        raw_sys_string = sys_string
                                                        break
                                                except:
                                                    continue
                                            else:
                                                sys_string = str(sub_system[0]).join(sys_string.split(str(sub_system)))
                
                                reduced_comps = sorted(flatten(eval(sys_string)))
                                central_ind = sys_comps[np.argmax(L_tot)]
                                positions = res['abspos'][reduced_comps]
                                central_pos = res['abspos'][central_ind]
                                separations = np.sqrt(np.sum((positions - central_pos)**2, axis=1))
                                non_zero_inds = np.where(separations>0)[0]
                                sep_list = sep_list + separations[non_zero_inds].tolist()
                                proximity_inds = np.argsort(separations)
                                p_ind = 2
                                L_list = L_list + [L_tot[proximity_inds[0]] + L_tot[proximity_inds[1]]]
                                while p_ind < len(proximity_inds):
                                    L_list.append(L_list[-1] + L_tot[proximity_inds[p_ind]])
                                    p_ind = p_ind + 1
                                    
                                invisible_sub_systems = np.where(sep_list < S_bins[bin_it-1])[0]
                                if len(invisible_sub_systems) > 0:
                                    separations_to_remove = np.array(sep_list)[invisible_sub_systems]
                                    for sep_rm in separations_to_remove:
                                        sep_list.remove(sep_rm)
                                    
                                    luminosity_to_remove = np.array(L_list)[invisible_sub_systems]
                                    for L_rm in luminosity_to_remove:
                                        L_list.remove(L_rm)
                        
                        if size == 1:
                            if str(sys_comps) not in All_unique_systems.keys():
                                try:
                                    sep_list = sep_list.tolist()
                                    All_unique_systems.update({str(sys_comps): sep_list})
                                    current_separations = current_separations + sep_list
                                except:
                                    All_unique_systems.update({str(sys_comps): sep_list})
                                    current_separations = current_separations + sep_list
                                #All_unique_systems_M.update({str(sys_comps): [M]})
                                if len(L_list)>0:
                                    All_unique_systems_L.update({str(sys_comps): [np.max(L_list)]})
                                    All_unique_systems_T.update({str(sys_comps): [time]})
                                else:
                                    All_unique_systems_L.update({str(sys_comps): [np.nan]})
                                    All_unique_systems_T.update({str(sys_comps): [np.nan]})
                            else:
                                try:
                                    All_unique_systems[str(sys_comps)].append(sep_list)
                                    current_separations = current_separations + sep_list
                                except:
                                    sep_list = sep_list.tolist()
                                    All_unique_systems[str(sys_comps)].append(sep_list)
                                    current_separations = current_separations + sep_list
                                #All_unique_systems_M[str(sys_comps)].append(M)
                                if len(L_list)>0:
                                    All_unique_systems_L[str(sys_comps)].append(np.max(L_list))
                                    All_unique_systems_T[str(sys_comps)].append(time)
                                else:
                                    All_unique_systems_L[str(sys_comps)].append(np.nan)
                                    All_unique_systems_T[str(sys_comps)].append(np.nan)
                        else:
                            if len(sep_list)>0:
                                sep_update.update({str(sys_comps): sep_list})
                            if len(L_list)>0:
                                L_update.update({str(sys_comps): [np.max(L_list)]})
                                T_update.update({str(sys_comps): [time]})
                            else:
                                L_update.update({str(sys_comps): [np.nan]})
                                T_update.update({str(sys_comps): [np.nan]})
                                
                        if len(sys_comps)>1:
                            L_tots.append(mean_L)
                    else:
                        res['n'][sys_comps] = 0
                        res['n'][multi_ind] = 0
        
                if size == 1:
                    L_tot_hist, bins = np.histogram(L_tots, bins=L_bins)
                    Luminosities.append(np.array(L_tot_hist))
                    if args.verbose_printing != 'False':
                        print("AFTER REDEFINING SYSTEMS WITH INVISIBLE COMPONENTS:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds])))
                        print("TOTAL NUMBER OF STARS =", str(missing_stars + np.sum(res['n'][top_inds])))
                    
                    ns = len(np.where(res['n'][top_inds]==1)[0])
                    nb = len(np.where(res['n'][top_inds]==2)[0])
                    nt = len(np.where(res['n'][top_inds]==3)[0])
                    nq = len(np.where(res['n'][top_inds]==4)[0])
                    nq5 = len(np.where(res['n'][top_inds]==5)[0])
                    ns6 = len(np.where(res['n'][top_inds]==6)[0])
                    ns7 = len(np.where(res['n'][top_inds]==7)[0])
                    n_systems.append([ns,nb,nt,nq,nq5,ns6,ns7])
                    cf = (nb+nt*2+nq*3+nq5*4+ns6*5+ns7*6)/(ns+nb+nt+nq+nq5+ns6+ns7)
                    if args.verbose_printing != 'False':
                        print("CF =", cf)
                    CF_per_bin.append(cf)
                    
            if size == 1:
                Separations.append(current_separations)
                N_sys_total.append(n_systems[0])
        try:
            if size != 1:
                L_tot_hist, bins = np.histogram(L_tots, bins=L_bins)
                if args.verbose_printing != 'False':
                    print("AFTER REDEFINING SYSTEMS WITH INVISIBLE COMPONENTS:" + str(len(np.where(res['n'][top_inds]==1)[0])) + ':' + str(len(np.where(res['n'][top_inds]==2)[0])) + ':' + str(len(np.where(res['n'][top_inds]==3)[0])) + ':' + str(len(np.where(res['n'][top_inds]==4)[0])) + ':' + str(len(np.where(res['n'][top_inds]==5)[0])) + ':' + str(len(np.where(res['n'][top_inds]==6)[0])) + ':' + str(len(np.where(res['n'][top_inds]==7)[0])) + '=' + str(np.sum(res['n'][top_inds])))
                    print("TOTAL NUMBER OF STARS =", str(missing_stars + np.sum(res['n'][top_inds])))

                if np.max(res['n']) > 6:
                    import pdb
                    pdb.set_trace()

                ns = len(np.where(res['n'][top_inds]==1)[0])
                nb = len(np.where(res['n'][top_inds]==2)[0])
                nt = len(np.where(res['n'][top_inds]==3)[0])
                nq = len(np.where(res['n'][top_inds]==4)[0])
                nq5 = len(np.where(res['n'][top_inds]==5)[0])
                ns6 = len(np.where(res['n'][top_inds]==6)[0])
                ns7 = len(np.where(res['n'][top_inds]==7)[0])
                n_systems.append([ns,nb,nt,nq,nq5,ns6,ns7])
                cf = (nb+nt*2+nq*3+nq5*4+ns6*5+ns7*6)/(ns+nb+nt+nq+nq5+ns6+ns7)
                if args.verbose_printing != 'False':
                    print("CF =", cf)
                CF_per_bin.append(cf)
                del ns
                del nb
                del nt
                del nq
                del nq5
                del ns6
                del ns7

                if rank > 0:
                    #print("Rank", rank, "All_unique_systems_L:", All_unique_systems_L)
                    data = {'current_separations': current_separations, 'CF_per_bin': CF_per_bin, 'n_systems': n_systems, 'L_tot_hist': L_tot_hist, 'All_unique_systems': sep_update, 'All_unique_systems_L': L_update, 'All_unique_systems_T': T_update}
                    CW.send(data, dest=0, tag=rank)
        except:
            pass
            
        sys.stdout.flush()
        CW.Barrier()
        
        if rank == 0:
            print("CF_per_bin before receiving data form other ranks:", CF_per_bin)
            #print("N_sys_total before receiving data form other ranks:", N_sys_total)
            
            for rit_recv in range(1,size):
                data = CW.recv(source=rit_recv, tag=rit_recv)
                Separations.append(data['current_separations'])
                CF_per_bin.append(data['CF_per_bin'][0])
                n_systems.append(data['n_systems'][0])
                Luminosities.append(np.array(data['L_tot_hist']))
                print("From rank:", rank, "received:", [data['All_unique_systems'], data['All_unique_systems_L'], data['All_unique_systems_T']])
                for key in data['All_unique_systems'].keys():
                    if key not in All_unique_systems.keys():
                        All_unique_systems.update({key: [data['All_unique_systems'][key]]})
                    else:
                        All_unique_systems[key].append(data['All_unique_systems'][key])
                for key in data['All_unique_systems_L'].keys():
                    if key not in All_unique_systems_L.keys():
                        All_unique_systems_L.update({key: [data['All_unique_systems_L'][key]]})
                        All_unique_systems_T.update({key: [data['All_unique_systems_T'][key]]})
                    else:
                        All_unique_systems_L[key].append(data['All_unique_systems_L'][key])
                        All_unique_systems_T[key].append(data['All_unique_systems_T'][key])
                else:
                    pass
            
            CF_Array_Full.append(CF_per_bin)
            N_sys_total.append(n_systems)
            #print("All_unique_systems:", All_unique_systems)
            
            #Write pickle
            pickle_file = savedir + args.pickled_file
            file = open(pickle_file, 'wb')
            pickle.dump((Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities),file)
            print('updated pickle', pickle_file, "for time_it", time_it, "of", end_time_ind+1)
            file.close()
            
        sys.stdout.flush()
        CW.Barrier()
    
#=====================================================
#Create plots below

if rank == 0:
    try:
        pickle_file = savedir + args.pickled_file
        file = open(pickle_file, 'rb')
    except:
        pickle_file = args.pickled_file
        file = open(pickle_file, 'rb')
    Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
    file.close()

    #CF calculated from all systems
    Summed_systems = np.sum(np.array(N_sys_total), axis=0)
    CF_top = Summed_systems[:,1] + Summed_systems[:,2]*2 + Summed_systems[:,3]*3 + Summed_systems[:,4]*4 + Summed_systems[:,5]*5 + Summed_systems[:,6]*6
    CF_bot = np.sum(Summed_systems, axis=1)
    CF_Total = CF_top/CF_bot

    plt.clf()
    try:
        plt.bar(bin_centers, CF_Total, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
        plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
    except:
        plt.bar(bin_centers[1:], CF_Total, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
        plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
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
        err_hist, err_bins = np.histogram(np.array(CF_Array_Full)[:,bit], bins=err_bins)
        width = err_bins[1] - err_bins[0]
        axes_dict[ax_label].bar(err_centers, err_hist, width=width)
        #plt.bar(err_centers, err_hist, width=width)
        
        xs = np.ones(np.shape(np.array(CF_Array_Full)[:,bit]))*bin_centers[bit]
        #plt.scatter(xs, np.array(CF_Array_Full)[:,bit], alpha=0.5, color='b')
        median = np.median(np.array(CF_Array_Full)[:,bit])
        axes_dict[ax_label].axvline(x=median, label='median', color='r', alpha=0.25)
        #plt.axvline(x=median, label='median', color='r')
        #plt.scatter(bin_centers[bit], median, color='r')
        mean = np.mean(np.array(CF_Array_Full)[:,bit])
        axes_dict[ax_label].axvline(x=mean, label='mean', color='orange', alpha=0.25)
        #plt.axvline(x=mean, label='mean', color='orange')
        #plt.scatter(bin_centers[bit], median, color='orange')
        std = np.std(np.array(CF_Array_Full)[:,bit], ddof=1)
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
        #plt.savefig(savedir+"error_dist_bin_"+str(bit)+".png", format='png', bbox_inches='tight')
        standard_deviation = [median-(mean-std), (mean+std)-median]
        #plt.plot([bin_centers[bit], bin_centers[bit]], [mean-std, mean+std], color='black')
        CF_median.append(median)
        CF_err.append(standard_deviation)

    plt.savefig(savedir+"error_dist.pdf", format='pdf', bbox_inches='tight')
    #plt.ylim(bottom=0.0)
    #plt.savefig(savedir+"error_dist.png")
    CF_err = np.array(CF_err)

    plt.clf()
    try:
        plt.bar(bin_centers, CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
    except:
        plt.bar(bin_centers[1:], CF_median, yerr=CF_err.T, edgecolor='k', label="CF Simulations", width=0.25, alpha=0.5)
    plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, edgecolor='black', alpha=0.5, label="Tobin et al")
    plt.legend(loc='best')
    plt.xlabel('Log Separation (AU)')
    plt.ylabel('Companion Frequency')
    plt.xlim([bin_centers[0]-0.25,bin_centers[-1]+0.25])
    plt.ylim(bottom=0.0)
    plt.savefig(savedir+'Median_companion_frequency.pdf', format='pdf', bbox_inches='tight')
    print('created Median_companion_frequency.pdf')

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

    #Iterate over systems and make histogram of masses of the high luminosity systems.
    """
    plt.clf()
    n_lines = 8
    c = np.arange(1, n_lines + 1)
    cmap = mpl.cm.get_cmap('jet', n_lines)
    dummie_cax = plt.scatter(c-100, c, c=c, cmap=cmap)
    for key in All_unique_systems.keys():
        if len(eval(key))>1:
            c_it = len(eval(key))
            if c_it > 8:
                c_it = 8
            x = np.arange(len(All_unique_systems[key]))
            plt.semilogy(x,All_unique_systems[key],c=cmap(c_it-1),alpha=0.5, lw=0.5)
    plt.colorbar(dummie_cax, ticks=c).set_label("Number of components")
    plt.xlabel("No. timesteps")
    plt.ylabel("Log Separation (AU)")
    plt.xlim(left=0.0)
    plt.savefig(savedir + args.figure_prefix + "Separation_unique_systems.jpg")

    n_cut_range = range(5,8+1)
    for n_cut in n_cut_range:
        Median_separation = []
        for key in All_unique_systems.keys():
            if len(eval(key))<n_cut:
                med_sep = np.median(All_unique_systems[key])
                Median_separation.append(med_sep)
        plt.clf()
        Median_hist, bins = np.histogram(Median_separation, bins=S_bins)
        plt.bar(((np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2), Median_hist, width=0.25, fill=False, edgecolor='black')
        plt.xlabel('Median Log Separation (AU)')
        plt.savefig(savedir + args.figure_prefix + "Median_Separation_N_less_than_"+str(n_cut)+"_Unique_Systems_hist.jpg")

    N_components = []
    for key in All_unique_systems_L.keys():
        components = eval(key)
        N_system = len(components)
        N_components.append(N_system)
    bins = np.arange(2,25)-0.5
    comp_hist, bins = np.histogram(N_components, bins=bins)
    plt.clf()
    plt.bar((bins[:-1]+bins[1:])/2, comp_hist, width=1, fill=False, edgecolor='black')
    plt.xlabel('Number of Components')
    plt.ylabel('Frequency')
    plt.savefig(savedir + args.figure_prefix + "Hist_N_components.jpg")
    """
    Mean_L = []
    for key in All_unique_systems_L.keys():
        Total_Luminosities = []
        non_nan_inds = np.where(np.isnan(All_unique_systems_T[key])==False)[0]
        if len(non_nan_inds) == 0 or len(eval(key))<2:
            life_time = 0
        else:
            life_time = (np.array(All_unique_systems_T[key])[non_nan_inds][-1] - np.array(All_unique_systems_T[key])[non_nan_inds][0])*units['time_unit'].in_units('yr')
        if life_time > args.lifetime_threshold:
            for L in All_unique_systems_L[key]:
                Total_Luminosities.append(np.max(L))
            if len(Total_Luminosities) > 0:
                Mean_L.append(np.mean(Total_Luminosities))
        #else:
        #    print("System", key, "is too short lived")
        
    
    L_tot_hist, bins = np.histogram(Mean_L, bins=L_bins)
    plt.clf()
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_tot_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Simulation")
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Tobin et al (2016)")
    plt.legend(loc='best')
    plt.xlabel('Mean Luminosty (log(L))')
    plt.ylabel('Number')
    #plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
    plt.ylim(bottom=0.0)
    plt.savefig(savedir + args.figure_prefix + 'Mean_luminosty_dist_of_unique_systems_with_L_phot.pdf', format='pdf', bbox_inches='tight')

    #plot max luminosities
    Max_L = []
    for key in All_unique_systems_L.keys():
        Total_Luminosities = []
        non_nan_inds = np.where(np.isnan(All_unique_systems_T[key])==False)[0]
        if len(non_nan_inds) == 0 or len(eval(key))<2:
            life_time = 0
        else:
            life_time = (np.array(All_unique_systems_T[key])[non_nan_inds][-1] - np.array(All_unique_systems_T[key])[non_nan_inds][0])*units['time_unit'].in_units('yr')
        if life_time > args.lifetime_threshold:
            for L in All_unique_systems_L[key]:
                Total_Luminosities.append(np.max(L))
            if len(Total_Luminosities) > 0:
                Max_L.append(np.max(Total_Luminosities))
        #else:
        #    print("System", key, "is too short lived")
    L_tot_hist, bins = np.histogram(Max_L, bins=L_bins)
    plt.clf()
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_tot_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Simulation")
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=(np.log10(L_bins[1])-np.log10(L_bins[0])), edgecolor='black', alpha=0.5, label="Tobin et al (2016)")
    plt.legend(loc='best')
    plt.xlabel('Max Luminosty (log(L))')
    plt.ylabel('Number')
    #plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
    plt.ylim(bottom=0.0)
    plt.savefig(savedir + args.figure_prefix + 'Max_luminosty_dist_of_unique_systems_with_L_phot.pdf', format='pdf', bbox_inches='tight')

    '''
    plt.clf()
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), np.sum(Luminosities, axis=0), width=0.31, edgecolor='black', label="Sum of luminosity hisotgrams")
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=0.31, edgecolor='black', fill=False, label="Tobin et al (2016)")
    plt.legend(loc='best')s
    plt.xlabel('Luminosty (log(L))')
    plt.ylabel('Number of instances')
    plt.yscale('log')
    plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
    plt.savefig(savedir+"Sum_of_Luminosity_histograms.jpg")
    '''
    '''
    plt.clf()
    used_sink_inds = []
    cmap = plt.get_cmap("tab10")
    window = 50
    for key in All_unique_systems.keys():
        append_bools = []
        for ind in eval(key):
            if ind not in flatten(used_sink_inds):
                append_bools.append(True)
            else:
                append_bools.append(False)
        if False not in append_bools:
            color = None
            used_sink_inds.append(eval(key))
        elif True not in append_bools:
            existing_ind = eval(key)[append_bools.index(False)]
            for sys_ind in used_sink_inds:
                if existing_ind in sys_ind:
                    color_ind = used_sink_inds.index(sys_ind)
                    break
            color = cmap(color_ind)
            print(key, "already in list")
        else:
            existing_ind = eval(key)[append_bools.index(False)]
            for sys_ind in used_sink_inds:
                if existing_ind in sys_ind:
                    true_inds = np.where(np.array(append_bools) == True)[0]
                    sys_ind.append(eval(key)[true_inds[0]])
                    color_ind = used_sink_inds.index(sys_ind)
                    break
            color = cmap(color_ind)
    
        system_times = np.array(All_unique_systems_T[key]).T[0]
        system_seps = np.array(All_unique_systems[key]).T[0]
        if np.max(system_seps) > 10000:
            remove_inds = np.argwhere(np.array(system_seps)>10000).T[0]
            system_seps[remove_inds[0]] = np.nan
            system_times[remove_inds[0]] = np.nan
        if len(system_times) == len(system_seps):
            lifetime = (np.array(All_unique_systems_T[key])[-1] - np.array(All_unique_systems_T[key])[0])*units['time_unit'].in_units('yr')
            print("lifetime for is", lifetime[0], "for", key)
            if lifetime > 0.0:
                moving_time = []
                moving_sep = []
                moving_upper_sep = []
                moving_lower_sep = []
                for smooth_ind in range(window, len(system_times)-window):
                    ave_time = np.mean(system_times[smooth_ind-window:smooth_ind+window])
                    ave_sep = np.mean(system_seps[smooth_ind-window:smooth_ind+window])
                    upper_sep = np.max(system_seps[smooth_ind-window:smooth_ind+window])
                    lower_sep = np.min(system_seps[smooth_ind-window:smooth_ind+window])
                    moving_time.append(ave_time)
                    moving_sep.append(ave_sep)
                    moving_upper_sep.append(upper_sep)
                    moving_lower_sep.append(lower_sep)
            if color == None:
                plt.semilogy(moving_time, moving_sep)
                plt.fill_between(moving_time, moving_lower_sep, moving_upper_sep, alpha=0.2)
            else:
                plt.semilogy(moving_time, moving_sep, color=color)
                plt.fill_between(moving_time, moving_lower_sep, moving_upper_sep, alpha=0.2, color=color)
    for bin_bound in S_bins:
        plt.axhline(y=bin_bound, color='k')
    plt.xlabel("Time ($t_{ff}$)")
    plt.ylabel("Log Separation (AU)")
    plt.savefig("Separation_evolution.jpg")
    '''
  
sys.stdout.flush()
CW.Barrier()
