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

window = yt.YTQuantity(args.time_intergration_window, 'yr')
luminosity_lower_limit = 0.01
accretion_limit = 1.e-7
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
    current_separations = []

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
    
    res = m.multipleAnalysis(S,nmax=6,cutoff=1e4, nmax=6, cyclic=False, Grho=Grho)
    time = S._time * units['time_unit'].in_units('yr')
    
    s = np.where((res['n']==1) & (res['topSystem']==True))[0]
    b = np.where((res['n']==2) & (res['topSystem']==True))[0]
    t = np.where((res['n']==3) & (res['topSystem']==True))[0]
    q = np.where((res['n']==4) & (res['topSystem']==True))[0]
    q5 = np.where((res['n']==5) & (res['topSystem']==True))[0]
    s6 = np.where((res['n']==6) & (res['topSystem']==True))[0]
    multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist())

    #CALCULATE Photospheric Luminosity
    radius = yt.YTQuantity(2.0, 'rsun')
    temperature = yt.YTQuantity(3000, 'K')
    L_phot = yt.units.stefan_boltzmann_constant * 4*np.pi*radius.in_units('cm')**2 * temperature**4
    
    #Find indexes of time window
    min_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time+window))
    max_time_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time-window))
    global_ind = np.argmin(abs(global_data['time'][:,0]*units['time_unit'].in_units('yr') - time))
    
    #Find Single stars that are within the limonsity limits
    dM = (global_data['m'][max_time_ind,s] - global_data['m'][min_time_ind,s])*units['mass_unit'].in_units('msun')
    dt = (global_data['time'][max_time_ind,s] - global_data['time'][min_time_ind,s])*units['time_unit'].in_units('yr')
    M_dot = dM/dt
    M = yt.YTArray(global_data['m'][global_ind,s]*units['mass_unit'].in_units('msun'), 'Msun')
    f_acc = 0.5
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    s_inds = np.where((L_tot>luminosity_lower_limit)&(M_dot>accretion_limit))[0] #&(L_tot<35.0)
    print("set res['n'] for invisible singles to 0")
    invisible_stars = list(set(s).symmetric_difference(s[s_inds]))
    res['n'][invisible_stars] = 0
    #Save single stars Luminosities
    
    #Filter out companions that are outside of the detection limits:
    Companion_inds = np.where((res['n']==1) & (res['topSystem']==False))[0]
    dM = (global_data['m'][max_time_ind,Companion_inds] - global_data['m'][min_time_ind,Companion_inds])*units['mass_unit'].in_units('msun')
    dt = (global_data['time'][max_time_ind,Companion_inds] - global_data['time'][min_time_ind,Companion_inds])*units['time_unit'].in_units('yr')
    M_dot = dM/dt
    M = yt.YTArray(global_data['m'][global_ind,Companion_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    f_acc = 0.5
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    detectable_companions = np.where((L_tot>luminosity_lower_limit))#&(M_dot>accretion_limit))[0]
    detectable_inds = Companion_inds[detectable_companions]
    
    for multi_ind in multi_inds:
        sys_comps = losi(multi_ind, res)
        sys_comps = sorted(flatten(sys_comps))
        detectable_components = list(set(sys_comps).intersection(detectable_inds))
        sys_comps = sorted(detectable_components)

        dM = (global_data['m'][max_time_ind,sys_comps] - global_data['m'][min_time_ind,sys_comps])*units['mass_unit'].in_units('msun')
        dt = (global_data['time'][max_time_ind,sys_comps] - global_data['time'][min_time_ind,sys_comps])*units['time_unit'].in_units('yr')
        M_dot = dM/dt
        M = yt.YTArray(global_data['m'][global_ind,sys_comps]*units['mass_unit'].in_units('msun'), 'Msun')
        f_acc = 0.5
        L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
        L_tot = L_acc.in_units('Lsun')
        mean_L = np.sum(L_tot)
        mean_M_dot = np.sum(M_dot)
        
        #if (mean_L+std_L) >0.04 and (mean_L-std_L) < 32.0 and mean_M_dot+std_M_dot>1.e-7:
        if mean_L > luminosity_lower_limit and mean_M_dot>accretion_limit: #&(L_tot<35.0)
            sys_comps = losi(multi_ind, res)
            sys_comps = sorted(flatten(sys_comps))
        
            #Check whether there is undetected companions
            n_detectable = len(detectable_components)
            res['n'][multi_ind] = n_detectable
            if n_detectable == len(sys_comps):
                #Check whether hierarchial. If it is it should be pretty easy for code:
                #res['n'][multi_ind] = n_detectable
                system_structure = losi(multi_ind, res)
            else:
                system_structure = losi(multi_ind, res)
                if n_detectable == 0:
                    #res['n'][multi_ind] = 0
                    print("none of the components are detected")
                elif n_detectable == 1:
                    sys_comps = detectable_components
                elif n_detectable == 2:
                    system_structure = detectable_components
                    sys_comps = sorted(flatten(system_structure))
                    
                    dM = (global_data['m'][max_time_ind,sys_comps] - global_data['m'][min_time_ind,sys_comps])*units['mass_unit'].in_units('msun')
                    dt = (global_data['time'][max_time_ind,sys_comps] - global_data['time'][min_time_ind,sys_comps])*units['time_unit'].in_units('yr')
                    M_dot = dM/dt
                    M = yt.YTArray(global_data['m'][global_ind,sys_comps]*units['mass_unit'].in_units('msun'), 'Msun')
                    f_acc = 0.5
                    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
                    L_tot = L_acc.in_units('Lsun')
                    mean_L = np.sum(L_tot)
                    mean_M_dot = np.sum(M_dot)
                    
                else:# len(sys_comps) - n_detectable == 1:
                    invisible_components = list(set(sys_comps).symmetric_difference(detectable_components))
                    sys_string = str(system_structure)
                    for inv_comp in invisible_components:
                        inv_string =str(inv_comp)
                        split_string = sys_string.split(inv_string)
                        sys_string = ''.join(split_string)
                        
                    reduced = False
                    while reduced == False:
                        if '[, ' in sys_string:
                            sys_string = ''.join(sys_string.split('[, '))[:-1*len(sys_string.split('[, '))+1]
                        if ' ],' in sys_string:
                            sys_string = ''.join(sys_string.split(' ],'))[1:]
                        if ', []' in sys_string:
                            sys_string = ''.join(sys_string.split(', []'))
                        if '[], ' in sys_string:
                            sys_string = ''.join(sys_string.split('[], '))
                        try:
                            system_structure = eval(sys_string)
                            reduced = True
                        except:
                            print("sys_string", sys_string, "is still not reduced")
                    
                    sys_comps = sorted(flatten(system_structure))
                    
                    dM = (global_data['m'][max_time_ind,sys_comps] - global_data['m'][min_time_ind,sys_comps])*units['mass_unit'].in_units('msun')
                    dt = (global_data['time'][max_time_ind,sys_comps] - global_data['time'][min_time_ind,sys_comps])*units['time_unit'].in_units('yr')
                    M_dot = dM/dt
                    M = yt.YTArray(global_data['m'][global_ind,sys_comps]*units['mass_unit'].in_units('msun'), 'Msun')
                    f_acc = 0.5
                    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
                    L_tot = L_acc.in_units('Lsun')
                    mean_L = np.sum(L_tot)
                    mean_M_dot = np.sum(M_dot)
                            
            #Now calculate separations accoridng to Tobin
            if n_detectable == 1 or n_detectable == 0:
                sys_comps = detectable_components
                sep_list = []
                L_list = []
            elif is_hierarchical(system_structure):
                central_ind = sys_comps[np.argmax(L_tot)]
                positions = res['abspos'][sys_comps]
                central_pos = res['abspos'][central_ind]
                separations = np.sqrt(np.sum((positions - central_pos)**2, axis=1))
                non_zero_inds = np.where(separations>0)[0]
                sep_list = separations[non_zero_inds].tolist()
                proximity_inds = np.argsort(separations)
                p_ind = 2
                L_list = [L_tot[proximity_inds[0]] + L_tot[proximity_inds[1]]]
                while p_ind < len(proximity_inds):
                    L_list.append(L_list[-1] + L_tot[proximity_inds[p_ind]])
                    p_ind = p_ind + 1
            else:
                sep_list = []
                L_list = []
                if len(flatten(system_structure)) == 4 and np.shape(system_structure) == (2,2):
                    primary_positions = []
                    for sys in system_structure:
                        binary_positions = res['abspos'][sys]
                        binary_sep = np.sqrt(np.sum((binary_positions[0] -  binary_positions[1])**2))
                        sep_list.append(binary_sep)
                        sys_comps_ind = []
                        for comp in sys:
                            sys_comps_ind.append(sys_comps.index(comp))
                        primary_ind = sys[np.argmax(np.array(L_tot)[sys_comps_ind])]
                        primary_positions.append(res['abspos'][primary_ind])
                        L_list.append(np.sum(L_tot[sys_comps_ind]))
                    quad_sep = np.sqrt(np.sum((primary_positions[0] -  primary_positions[1])**2))
                    sep_list.append(quad_sep)
                    L_list.append(np.sum(L_tot))
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
                                    binary_sep = np.sqrt(np.sum((binary_positions[0] -  binary_positions[1])**2))
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
            
            if str(sys_comps) not in All_unique_systems.keys():
                try:
                    sep_list = sep_list.tolist()
                    All_unique_systems.update({str(sys_comps): [sep_list]})
                    current_separations = current_separations + sep_list
                except:
                    All_unique_systems.update({str(sys_comps): [sep_list]})
                    current_separations = current_separations + sep_list
                All_unique_systems_M.update({str(sys_comps): [M]})
                All_unique_systems_L.update({str(sys_comps): [L_list]})
            else:
                try:
                    All_unique_systems[str(sys_comps)].append(sep_list)
                    current_separations = current_separations + sep_list
                except:
                    sep_list = sep_list.tolist()
                    All_unique_systems[str(sys_comps)].append(sep_list)
                    current_separations = current_separations + sep_list
                All_unique_systems_M[str(sys_comps)].append(M)
                All_unique_systems_L[str(sys_comps)].append(L_list)
            L_tots.append(mean_L)
    
    #Plot Luminosity histogram
    L_tot_hist, bins = np.histogram(L_tots, bins=L_bins)
    Luminosities.append(np.array(L_tot_hist))
    '''
    plt.clf()
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_tot_hist, width=0.31, edgecolor='black', label="Simulation")
    plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=0.31, edgecolor='black', fill=False, label="Tobin et al (2016)")
    plt.legend(loc='best')
    plt.xlabel('Luminosty (log(L))')
    plt.ylabel('Number')
    plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
    plt.text(-1.4, 20.5, 't='+str(int(np.round(time)))+'$yr$')
    plt.savefig(savedir + args.figure_prefix + 'luminosty_dist_' + str(int(np.round(time))) + '.jpg')
    print('created luminosty_dist_' + str(int(np.round(time))) + '.jpg')
    '''
    
    #Now bin data and create plots:
    cf_array = []
    n_systems = []
    print("CF CALCULATION PER BIN IS WRONG, IGNORE FOR NOW.")
    s = np.where(res['n']==1)[0]
    for bin_it in range(1,len(S_bins)):
        bin_inds = np.where((current_separations>=S_bins[bin_it-1])&(current_separations<S_bins[bin_it]))[0]
        CF_single_inds = np.where((current_separations<S_bins[bin_it-1])|(current_separations>=S_bins[bin_it]))[0]
        
        b = np.where(res['n'][bin_inds]==2)[0]
        b_single = np.where(res['n'][CF_single_inds]==2)[0]
        nb = np.count_nonzero(res['topSystem'][bin_inds][b])
        
        t = np.where(res['n'][bin_inds]==3)[0]
        t_single = np.where(res['n'][CF_single_inds]==3)[0]
        nt = np.count_nonzero(res['topSystem'][bin_inds][t])
        
        q = np.where(res['n'][bin_inds]==4)[0]
        q_single = np.where(res['n'][CF_single_inds]==4)[0]
        nq = np.count_nonzero(res['topSystem'][bin_inds][q])
        
        q5 = np.where(res['n'][bin_inds]==5)[0]
        q5_single = np.where(res['n'][CF_single_inds]==5)[0]
        nq5 = np.count_nonzero(res['topSystem'][bin_inds][q5])
        
        s6 = np.where(res['n'][bin_inds]==6)[0]
        s6_single = np.where(res['n'][CF_single_inds]==6)[0]
        ns6 = np.count_nonzero(res['topSystem'][bin_inds][s6])
        
        ns = np.count_nonzero(res['topSystem'][s]) + len(b_single) + len(t_single) + len(q_single) + len(q5_single) + len(s6_single)
        cf = (nb+nt*2+nq*3+nq5*4+ns6*5)/(ns+nb+nt+nq+nq5+ns6)
        n_systems.append([ns,nb,nt,nq,nq5,ns6])
        cf_array.append(cf)
        
    Separations.append(current_separations)
    Times.append(time)
    
    #Create histograms:
    '''
    plt.clf()
    Sep_hist, bins = np.histogram(current_separations, bins=S_bins)
    plt.bar(((np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2), Sep_hist, width=0.25, fill=False, edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.xlim([1,4])
    plt.savefig(savedir + args.figure_prefix + 'separation_histogram_' + str(int(np.round(time))) + '.jpg')
    print('created separation_histogram_' + str(int(np.round(time))) + '.jpg')
    
    plt.clf()
    plt.bar(((np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2), cf_array, width=0.25, fill=False, edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel('Companion Frequency')
    plt.xlim([1,4])
    plt.ylim([0.0, 0.25])
    plt.text(1.1, 0.23, 't='+str(int(np.round(time)))+'$yr$')
    plt.savefig(savedir + args.figure_prefix + 'companion_frequency_' + str(int(np.round(time))) + '.jpg')
    print('created companion_frequency_' + str(int(np.round(time))) + '.jpg')
    '''
    CF_Array_Full.append(cf_array)
    N_sys_total.append(n_systems)
    print("Calculated CFs of:", cf_array)
    
    #Update saved arrays:
   
#Print Calculate summed values
N_sys_total = np.array(N_sys_total)
N_total = np.sum(N_sys_total,axis=0)

N_sys_sum = np.sum(N_sys_total, axis=0)
CF_bottom = np.sum(N_sys_sum,axis=1)
CF_top = N_sys_sum[:, 1] + N_sys_sum[:, 2]*2 + N_sys_sum[:, 3]*3 + N_sys_sum[:, 4]*4 + N_sys_sum[:, 5]*5

#CF_bottom = np.sum(N_total, axis=1)
#CF_top = N_sys_sum[:, 1] + N_sys_sum[:, 2]*2 + N_sys_sum[:, 3]*3 + N_sys_sum[:, 4]*4 + N_sys_sum[:, 5]*5 + N_sys_sum[:, 6]*6 + N_sys_sum[:, 7]*7 + N_sys_sum[:, 8]*8
CF_total = CF_top/CF_bottom

plt.clf()
plt.bar(((np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2), CF_total, width=0.25, fill=False, edgecolor='black')
plt.xlabel('Separation')
plt.ylabel('Companion Frequency')
plt.xlim([1,4])
#plt.ylim([0.0, 0.25])
plt.savefig(savedir +  args.figure_prefix + 'Total_companion_frequency_.jpg')
print('created Total_companion_frequency.jpg')

CF_Array_Full = np.array(CF_Array_Full)
CF_sum = np.sum(CF_Array_Full, axis=0)

plt.clf()
plt.bar(((np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2), CF_sum, width=0.25, fill=False, edgecolor='black')
plt.xlabel('Separation')
plt.ylabel('Companion Frequency')
plt.xlim([1,4])
plt.savefig(savedir + args.figure_prefix + 'Sum_companion_frequency_.jpg')
print('created Sum_companion_frequency.jpg')

#Iterate over systems and make histogram of masses of the high luminosity systems.

All_masses = np.array([])
N_components = []
for key in All_unique_systems_L.keys():
    if np.log10(np.mean(All_unique_systems_L[key])) > 0.7:
        try:
            All_masses = np.append(All_masses, All_unique_systems_M[key][-1].value)
            N_components.append(len(All_unique_systems_M[key][-1]))
        except:
            All_masses = np.append(All_masses, All_unique_systems_M[key])
            N_components.append(len(All_unique_systems_M[key]))
            
comp_bins = np.linspace(1.5,5.5,5)
Comp_hist, bins = np.histogram(N_components, bins=comp_bins)
plt.clf()
plt.bar(((comp_bins[:-1]+comp_bins[1:])/2), Comp_hist, width=1.0, edgecolor='black', alpha=0.5, label="Simulation")
plt.xlabel('Number of Components')
plt.ylabel('Number')
plt.xlim([comp_bins[0],comp_bins[-1]])
plt.savefig(savedir + args.figure_prefix + 'N_components_of_system_mean_L_greater_than_0.7.jpg')

            
mass_bins = np.logspace(-2,0.5,11)
Mass_hist, bins = np.histogram(All_masses, bins=mass_bins)
plt.clf()
plt.bar(((np.log10(mass_bins[:-1])+np.log10(mass_bins[1:]))/2), Mass_hist, width=0.25, edgecolor='black', alpha=0.5, label="Simulation")
plt.xlabel('Final Mass (log(M))')
plt.ylabel('Number')
plt.xlim([np.log10(mass_bins[0]),np.log10(mass_bins[-1])])
plt.savefig(savedir + args.figure_prefix + 'Mass_histograms_of_system_mean_L_greater_than_0.7.jpg')

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
plt.ylabel("Separation (AU)")
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
    plt.xlabel('Median Separation (AU)')
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

Mean_L = []
for key in All_unique_systems_L.keys():
    if len(eval(key))>1:
        Mean_L.append(np.mean(np.array(All_unique_systems_L[key])))
L_tot_hist, bins = np.histogram(Mean_L, bins=L_bins)
plt.clf()
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_tot_hist, width=0.31, edgecolor='black', alpha=0.5, label="Simulation")
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=0.31, edgecolor='black', alpha=0.5, label="Tobin et al (2016)")
plt.legend(loc='best')
plt.xlabel('Mean Luminosty (log(L))')
plt.ylabel('Number')
plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
plt.savefig(savedir + args.figure_prefix + 'Mean_luminosty_dist_of_unique_systems_with_L_phot.jpg')

#plot max luminosities
Max_L = []
for key in All_unique_systems_L.keys():
    Max_L.append(np.mean(np.array(All_unique_systems_L[key])))
L_tot_hist, bins = np.histogram(Max_L, bins=L_bins)
plt.clf()
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), L_tot_hist, width=0.31, edgecolor='black', label="Simulation")
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=0.31, edgecolor='black', fill=False, label="Tobin et al (2016)")
plt.legend(loc='best')
plt.xlabel('Max Luminosty (log(L))')
plt.ylabel('Number')
plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
plt.savefig(savedir + args.figure_prefix + 'Max_luminosty_dist_of_unique_systems_with_L_phot.jpg')

plt.clf()
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), np.sum(Luminosities, axis=0), width=0.31, edgecolor='black', label="Sum of luminosity hisotgrams")
plt.bar(((np.log10(L_bins[:-1])+np.log10(L_bins[1:]))/2), Tobin_hist, width=0.31, edgecolor='black', fill=False, label="Tobin et al (2016)")
plt.legend(loc='best')
plt.xlabel('Luminosty (log(L))')
plt.ylabel('Number of instances')
plt.yscale('log')
plt.xlim([np.log10(L_bins[0]),np.log10(L_bins[-1])])
plt.savefig(savedir+"Sum_of_Luminosity_histograms.jpg")
