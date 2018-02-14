import binary_orbit as bo
import numpy as np
from astropy.time import Time
import csv
import matplotlib.pyplot as plt
from tqdm import tqdm

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", "--input_file", help="what file do you want to read in?", type=str)
    parser.add_argument("-bayes_f", "--bayes_file", help="where do you want to read out the calculated bayes factors?", type=str, default='bayes_factors.csv')
    parser.add_argument("-calc", "--calculate_bayes", type=str, default='False')
    parser.add_argument("-suffix", "--image_suffix", help="do you want to add a suffix to the saved images?", type=str, default='')
    args = parser.parse_args()
    return args

args = parse_inputs()

n_orb=int(1e4)
n_systems = int(1e2)
plotit=False
q_min = 0.05
my_orb = bo.random_orbits(n_orb=n_orb)
US_group_vel = 10.
UCL_group_vel = 4.
US_group_std = 1.3 * 3
UCL_group_std = 1.3 * 3
standard_std = {'F':1.08, 'G':0.67, 'K':1.43, 'M':2.27}# 2.0
astrophysical_std = 1.0 #Astrophysical radial velocity uncertainty

Object = []
Region = []
Coords = []
Membership_probability = []
IR_excess = []
Disk_visual = []
Magnitudes = []
No_obs = []
SB1_flag = []
SB2_flag = []
RV_variation = []
Temp_sptype = []
Pref_template = []
Obs_info = []
SNR = []
all_bayes = [[],[]]

RV_standard_info = {}

#Read in RV standard list
header = 0
with open('/Users/rajikak/Observational_Data/RV_standard_list.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            RV_standard_info[row[0]] = (float(row[5]), float(row[6]), float(row[7]))
        else:
            header = 1

print "Reading in current spreadsheet"
header = 0
with open(args.input_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if 'U4' in row[0]:
                row[0] = 'UCAC4' + row[0].split('U4')[-1]
            Object.append(row[0])
            Region.append(row[1])
            Coords.append((row[2], row[3]))
            Membership_probability.append(int(row[4]))
            IR_excess.append(row[5])
            Disk_visual.append(row[6])
            Magnitudes.append([float(row[7]), float(row[8]), float(row[9])])
            No_obs.append(int(row[10]))
            if row[11] == 'FALSE' or row[11] == 'False':
                SB1_flag.append(False)
            else:
                SB1_flag.append(True)
            if row[12] == 'FALSE' or row[12] == 'False':
                SB2_flag.append(False)
            else:
                SB2_flag.append(True)
            RV_variation.append(row[13])
            Pref_template.append(row[14])
            Temp_sptype.append(row[15])
            if len(row) > 16:
                Obs = np.array(row[16:])
                Obs = np.delete(Obs, np.where(Obs==''))
                #if len(Obs) > 5:
                #    Obs = np.reshape(Obs, (len(Obs)/5, 5))
                Obs = np.reshape(Obs, (len(Obs)/6, 6))
                for ind_obs in Obs:
                    if '/' in ind_obs[0]:
                        new_format = '20' + ind_obs[0].split('/')[-1] + '-' + ind_obs[0].split('/')[-2] + '-' + ("%02d" % int(ind_obs[0].split('/')[-3]))
                        ind_obs[0] = new_format
            else:
                Obs = np.array([])
            Obs_info.append(Obs)
        if header == 0:
            header = 1

Obj_bayes = np.nan*np.zeros(len(Object))
Obj_bayes_SpT = [[],[]]

#Read in currently calculated Bayes Factors:
header = 0
with open(args.bayes_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            ind = Object.index(row[0])
            Obj_bayes[ind] = float(row[2])
            if row[1] == 'US':
                all_bayes[0].append(float(row[2]))
                Obj_bayes_SpT[0].append(row[3])
            else:
                all_bayes[1].append(float(row[2]))
                Obj_bayes_SpT[1].append(row[3])
        else:
            header = 1

if args.calculate_bayes == 'False':
    all_bayes[0] = np.nan_to_num(np.array(all_bayes[0]))
    all_bayes[1] = np.nan_to_num(np.array(all_bayes[1]))
    all_bayes[0] = np.minimum(all_bayes[0],1e6)
    all_bayes[1] = np.minimum(all_bayes[1],1e6)

    ####### CALCULATED BAYES FACTORS HAVE BEEN READ IN

    gamma = np.linspace(0.0, 1.0, 1000)
    P_gamma = 1./(gamma*(1-gamma))
    gamma_prior = 1.0 #A constant for now.
    p_gamma_US = 1.0
    p_gamma_UCL = 1.0

    #US
    for bayes in all_bayes[0]:
        p_gamma_US *= gamma*bayes + (1.0-gamma)
    p_gamma_US = P_gamma*p_gamma_US
    p_gamma_US_norm = p_gamma_US/(np.trapz(p_gamma_US, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_US =  np.trapz(p_gamma_US_norm*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_US = np.trapz(p_gamma_US_norm*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_US = np.sqrt(mean_squared_US - mean_US*mean_US)

    #Calculate percentiles
    prev_delta = 1.
    for x_gamma in range(len(gamma)-2)[1:]:
        area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
        delta = np.abs(area - 0.16)
        if delta > prev_delta:
            US_16th = gamma[x_gamma-1]
            break
        prev_delta = delta
    prev_delta = 1.
    for x_gamma in range(len(gamma)-2)[1:]:
        area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
        delta = np.abs(area - 0.84)
        if delta > prev_delta:
            US_84th = gamma[x_gamma-1]
            break
        prev_delta = delta
    US_pos_error = US_84th - mean_US
    US_neg_error = mean_US - US_16th

    plt.clf()
    plt.plot(gamma, p_gamma_US_norm)
    #plt.title('US:' + "{0:5.1f}".format(mean_US*100.)+'$^{' + "{0:5.1f}".format(US_pos_error*100.)+'}_{' + "{0:5.1f}".format(US_neg_error*100.)+'}$%')
    plt.savefig('Binary_fraction_US'+args.image_suffix+'.png')

    #UCL
    for bayes in all_bayes[1]:
        p_gamma_UCL *= gamma*bayes + (1.0-gamma)
    p_gamma_UCL = P_gamma*p_gamma_UCL
    p_gamma_UCL_norm = p_gamma_UCL/(np.trapz(p_gamma_UCL, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_UCL =  np.trapz(p_gamma_UCL_norm*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_UCL = np.trapz(p_gamma_UCL_norm*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_UCL = np.sqrt(mean_squared_UCL - mean_UCL*mean_UCL)

    #Calculate percentiles
    prev_delta = 1.
    for x_gamma in range(len(gamma)-2)[1:]:
        area = np.trapz(p_gamma_UCL_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
        delta = np.abs(area - 0.16)
        if delta > prev_delta:
            UCL_16th = gamma[x_gamma-1]
            break
        prev_delta = delta
    prev_delta = 1.
    for x_gamma in range(len(gamma)-2)[1:]:
        area = np.trapz(p_gamma_UCL_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
        delta = np.abs(area - 0.84)
        if delta > prev_delta:
            UCL_84th = gamma[x_gamma-1]
            break
        prev_delta = delta
    UCL_pos_error = UCL_84th - mean_UCL
    UCL_neg_error = mean_UCL - UCL_16th

    plt.clf()
    plt.plot(gamma, p_gamma_UCL_norm)
    #plt.title('UCL:' + "{0:5.1f}".format(mean_UCL*100.)+'$^{' + "{0:5.1f}".format(UCL_pos_error*100.)+'}_{' + "{0:5.1f}".format(UCL_neg_error*100.)+'}$%')
    plt.savefig('Binary_fraction_UCL'+args.image_suffix+'.png')

    plt.clf()
    plt.plot(gamma, p_gamma_US_norm, label='Upper Scorpius')
    plt.plot(gamma, p_gamma_UCL_norm, label='Upper Centaurus-Lupus')
    plt.axvline(x=0.12)
    #plt.fill_between([11,13], np.zeros(len(gamma)), np.ones(len(gamma)), where=y2 >= y1, facecolor='green', interpolate=True, alpha=0.5)
    plt.legend(loc='best')
    plt.xlabel('Binary Fraction')
    plt.ylabel('Likelihood')

    #plt.title('US:' + "{0:5.1f}".format(mean_US*100.)+'$^{' + "{0:5.1f}".format(US_pos_error*100.)+'}_{' + "{0:5.1f}".format(US_neg_error*100.)+'}$%, UCL:' + "{0:5.1f}".format(mean_UCL*100.)+'$^{' + "{0:5.1f}".format(UCL_pos_error*100.)+'}_{' + "{0:5.1f}".format(UCL_neg_error*100.)+'}$%')
    plt.savefig('Composite_fractions' + args.image_suffix + '.png')

    ## per spectral type
    US_per_SpT = []
    UCL_per_SpT = []
    US_per_SpT_err = []
    UCL_per_SpT_err = []
    gamma_prior = 1.0 #A constant for now.
    p_gamma_US_F = 1.0
    p_gamma_UCL_F = 1.0
    p_gamma_US_G = 1.0
    p_gamma_UCL_G = 1.0
    p_gamma_US_K = 1.0
    p_gamma_UCL_K = 1.0
    p_gamma_US_M = 1.0
    p_gamma_UCL_M = 1.0

    for bayes in range(len(all_bayes[0])):
        if Obj_bayes_SpT[0][bayes] == 'F':
            p_gamma_US_F *= gamma*all_bayes[0][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[0][bayes] == 'G':
            p_gamma_US_G *= gamma*all_bayes[0][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[0][bayes] == 'K':
            p_gamma_US_K *= gamma*all_bayes[0][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[0][bayes] == 'M':
            p_gamma_US_M *= gamma*all_bayes[0][bayes] + (1.0-gamma)
    p_gamma_US_norm_F = p_gamma_US_F/(np.trapz(p_gamma_US_F, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_US_F =  np.trapz(p_gamma_US_norm_F*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_US_F = np.trapz(p_gamma_US_norm_F*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_US_F = np.sqrt(mean_squared_US_F - mean_US_F*mean_US_F)

    p_gamma_US_norm_G = p_gamma_US_G/(np.trapz(p_gamma_US_G, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_US_G =  np.trapz(p_gamma_US_norm_G*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_US_G = np.trapz(p_gamma_US_norm_G*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_US_G = np.sqrt(mean_squared_US_G - mean_US_G*mean_US_G)

    p_gamma_US_norm_K = p_gamma_US_K/(np.trapz(p_gamma_US_K, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_US_K =  np.trapz(p_gamma_US_norm_K*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_US_K = np.trapz(p_gamma_US_norm_K*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_US_K = np.sqrt(mean_squared_US_K - mean_US_K*mean_US_K)

    p_gamma_US_norm_M = p_gamma_US_M/(np.trapz(p_gamma_US_M, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_US_M =  np.trapz(p_gamma_US_norm_M*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_US_M = np.trapz(p_gamma_US_norm_M*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_US_M = np.sqrt(mean_squared_US_M - mean_US_M*mean_US_M)

    US_per_SpT.append(mean_squared_US_F)
    US_per_SpT_err.append(sigma_US_F)
    US_per_SpT.append(mean_squared_US_G)
    US_per_SpT_err.append(sigma_US_G)
    US_per_SpT.append(mean_squared_US_K)
    US_per_SpT_err.append(sigma_US_K)
    US_per_SpT.append(mean_squared_US_M)
    US_per_SpT_err.append(sigma_US_M)

    for bayes in range(len(all_bayes[1])):
        if Obj_bayes_SpT[1][bayes] == 'F':
            p_gamma_UCL_F *= gamma*all_bayes[1][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[1][bayes] == 'G':
            p_gamma_UCL_G *= gamma*all_bayes[1][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[1][bayes] == 'K':
            p_gamma_UCL_K *= gamma*all_bayes[1][bayes] + (1.0-gamma)
        if Obj_bayes_SpT[1][bayes] == 'M':
            p_gamma_UCL_M *= gamma*all_bayes[1][bayes] + (1.0-gamma)
    '''
    p_gamma_UCL_norm_F = p_gamma_UCL_F/(np.trapz(p_gamma_UCL_F, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_UCL_F =  np.trapz(p_gamma_UCL_norm_F*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_UCL_F = np.trapz(p_gamma_UCL_norm_F*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_UCL_F = np.sqrt(mean_squared_UCL_F - mean_UCL_F*mean_UCL_F)
    '''
    
    p_gamma_UCL_norm_G = p_gamma_UCL_G/(np.trapz(p_gamma_UCL_G, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_UCL_G =  np.trapz(p_gamma_UCL_norm_G*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_UCL_G = np.trapz(p_gamma_UCL_norm_G*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_UCL_G = np.sqrt(mean_squared_UCL_G - mean_UCL_G*mean_UCL_G)
    
    p_gamma_UCL_norm_K = p_gamma_UCL_K/(np.trapz(p_gamma_UCL_K, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_UCL_K =  np.trapz(p_gamma_UCL_norm_K*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_UCL_K = np.trapz(p_gamma_UCL_norm_K*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_UCL_K = np.sqrt(mean_squared_UCL_K - mean_UCL_K*mean_UCL_K)
    
    p_gamma_UCL_norm_M = p_gamma_UCL_M/(np.trapz(p_gamma_UCL_M, x=gamma, dx=(gamma[1]-gamma[0])))
    mean_UCL_M =  np.trapz(p_gamma_UCL_norm_M*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    mean_squared_UCL_M = np.trapz(p_gamma_UCL_norm_M*gamma*gamma, x=gamma, dx=(gamma[1]-gamma[0]))
    sigma_UCL_M = np.sqrt(mean_squared_UCL_M - mean_UCL_M*mean_UCL_M)
    
    UCL_per_SpT.append(np.nan)
    UCL_per_SpT_err.append(np.nan)
    UCL_per_SpT.append(mean_squared_UCL_G)
    UCL_per_SpT_err.append(sigma_UCL_G)
    UCL_per_SpT.append(mean_squared_UCL_K)
    UCL_per_SpT_err.append(sigma_UCL_K)
    UCL_per_SpT.append(mean_squared_UCL_M)
    UCL_per_SpT_err.append(sigma_UCL_M)

    x_rag = [0, 1.75, 2.33, 2.66, 4.33]
    x_rag_err = [1.0, 0.3, 1.0, 0.6, 0.5]
    y_rag_err = [0.3, 0.05, 0.025, 0.03, 0.05]
    y_rag = [0.75, 0.5, 0.45, 0.40, 0.35]
    x = [0,1,2,3,4]
    my_xticks = ['A','F','G','K','M']
    plt.clf()
    plt.xticks(x, my_xticks)
    plt.errorbar(x_rag, y_rag, xerr=x_rag_err, yerr=y_rag_err, fmt='+', color='k', label='Raghavan et al. (2010)')
    plt.errorbar(x, US_per_SpT, yerr=US_per_SpT_err, fmt='o', color='b', label='Upper Scorpius')
    plt.errorbar(x, UCL_per_SpT, yerr=UCL_per_SpT_err, fmt='o', color='r', label='Upper Centaurus-Lupus')
    plt.legend(loc='best')
    plt.xlabel('Spectral Type of template')
    plt.ylabel('Binary Fraction')
    plt.savefig('Per_SpT' + args.image_suffix + '.eps')


f = open(args.bayes_file, 'w')
f.write('Object,Region,Bayes_factor,SpT\n')
for obj_ind in range(len(Object)):
    if np.isnan(Obj_bayes[obj_ind]) == False:
        write_string = Object[obj_ind] + ',' + Region[obj_ind] + ',' + str(Obj_bayes[obj_ind]) + ',' + Temp_sptype[obj_ind][0:2] + '\n'
        f.write(write_string)

f.close()


if args.calculate_bayes == 'True':
    for obj in range(len(Object)):
        #Get M_1 from the template spectral type
        Pref_template_name = Pref_template[obj].split('_')[0]
        if Pref_template_name != '' and 'Y' in IR_excess[obj] and np.isnan(Obj_bayes[obj]):
            
            likelihoods = []
            single_likelihoods = []
        
            #Produces masses within +/- 10% of the mass of the template.
            #!!! Mike suggests a single mass.
            M_1 = (np.random.random(n_systems)*(RV_standard_info[Pref_template_name][1]-RV_standard_info[Pref_template_name][0])) + RV_standard_info[Pref_template_name][0]
            
            #Generates mass ratios with minium mass ratio of q_min (default 0.01?, should this be dependant on the primary mass? Because sometimes low mass ratios could give very low mass companions i.e. BD mass...)
            #!!! Mike suggests 0.05 due to brown dwarf desert.
            q = (np.random.random(n_systems)*(1-q_min)) + q_min
            
            #from Primary masses and mass ratios, secondary masses can get calculated
            M_2 = M_1 * q
            
            #Get dates of the observations of the object
            jds = Obs_info[obj][:,1].astype(np.float)
            jds_more = np.linspace(jds[0],jds[-1],1000)
            
            #get observed data, and add in the error in the standards in quadrature.
            #This relates to the spectrograph stability
            #There is also an astrophysical error due to these objects being rapid rotators etc.
            RV_standard_err = standard_std[Temp_sptype[obj][0]]
            err = np.sqrt(Obs_info[obj][:,3].astype(float)**2. + RV_standard_err**2. + astrophysical_std**2.)
            observed_rv = Obs_info[obj][:,2].astype(float)
            
            #IN A LOOP iterate over random orbits:
            for orb in tqdm(range(n_orb)):
                #FIXME: Figure out which velocity to use!
                if Region[obj] == 'US':
                    #v_group = np.random.normal(np.mean(observed_rv), np.sqrt(US_group_std**2 + RV_standard_err**2), n_systems)
                    v_group = np.random.normal(US_group_vel, np.sqrt(US_group_std**2 + RV_standard_err**2), n_systems)
                else:
                    #v_group = np.random.normal(np.mean(observed_rv), np.sqrt(UCL_group_std**2 + RV_standard_err**2), n_systems)
                    v_group = np.random.normal(UCL_group_vel, np.sqrt(UCL_group_std**2 + RV_standard_err**2), n_systems)

                #generate orbit?
                #!!! Find just one set of orbital parameters at at a time, and
                #scale the RVS. OR if you really want you can compute a, i etc
                #yourself and plug these into my_orb, but some RV scalign is still needed.
                rho, theta, normalised_vr = bo.binary_orbit(my_orb, jds, plot_orbit_no=orb)
                more_rho, more_theta, more_normalised_vr = bo.binary_orbit(my_orb, jds_more, plot_orbit_no=orb)
                for system in range(n_systems):
                    actual_vr = bo.scale_rv(normalised_vr, my_orb['P'][orb], M_1[system], M_2[system], my_orb['i'][orb], group_velocity=v_group[system])
                    actual_vr_more = bo.scale_rv(more_normalised_vr, my_orb['P'][orb], M_1[system], M_2[system], my_orb['i'][orb], group_velocity=v_group[system])
                    
                    this_likelihood = bo.calc_likelihood(actual_vr, observed_rv, err)
                    likelihoods.append(this_likelihood)
                    #THEN CALCULATE PROBABILITY OF BEING A SINGLE STAR
                    single_likelihoods.append(bo.calc_likelihood(v_group[system], observed_rv, err))
                    
                    if (this_likelihood > 0.95) and plotit:
                        #if np.sum(np.abs(observed_rv-actual_vr) < err) > 3.:
                        #plot orbit
                        print "PLOTTING ORBIT for maximum_likelihood =", this_likelihood
                        plt.clf()
                        plt.plot(jds_more, actual_vr_more)
                        plt.scatter(jds, actual_vr, color='k')
                        plt.errorbar(jds, observed_rv, yerr=err, fmt='o', color='r')
                        plt.show()
                        print "MOVING ONTO NEXT ORBIT"
            
            #THEN CALCULATE BAYES FACTOR
            bayes_factor = np.mean(likelihoods)/np.mean(single_likelihoods)
            print("Bayes Factor: {0:5.2f} for ".format(bayes_factor) + Object[obj])
            if np.isinf(bayes_factor):
                bayes_factor = np.nan_to_num(bayes_factor)
            if Region[obj] == 'US':
                all_bayes[0].append(bayes_factor)
            else:
                all_bayes[1].append(bayes_factor)
            Obj_bayes[obj] = bayes_factor
            f = open(args.bayes_file, 'a')
            write_string = Object[obj] + ',' + Region[obj] + ',' + str(bayes_factor) + '\n'
            f.write(write_string)
            f.close()