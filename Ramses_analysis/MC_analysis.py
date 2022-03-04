import glob
import numpy as np
import pickle

def calc_likelihood(scaled_rvs, observed_rvs, err):
    """
        This function calculates the likelihood of an orbit fitting the observed data
    """
    chi_squared = np.sum((observed_rvs - scaled_rvs)**2./err**2.)
    maximum_likelihood = np.exp(-chi_squared/2.)
    return maximum_likelihood # chi_squared

def chi_squared(expected, observed):
    chi_squared = np.sum(((observed - expected)**2)/expected)
    return chi_squared

pickle_files = sorted(glob.glob("*/*.pkl"))

#Read Tobin data:
file = open("Tobin_CF.pkl", 'rb')
S_bins, CF_per_bin_Tobin = pickle.load(file)
file.close()

CF_per_bin_Tobin = np.array(CF_per_bin_Tobin)

print("Read Tobin Data")

bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

#Read unbound absolute separation distribution:
abs_pickle = "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/Limits_on_L_and_Acc/Unbound/Using_Global_Pickle/para_test_both_limits_unbound.pkl"
file = open(abs_pickle, 'rb')
Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
file.close()

Abs_median = np.median(np.array(CF_Array_Full), axis=0)
mean = np.mean(np.array(CF_Array_Full), axis=0)
std = np.std(np.array(CF_Array_Full), axis=0)
Abs_err = np.array([Abs_median-(mean-std), (mean+std)-Abs_median])
Abs_median = np.array(Abs_median)

print("Calculated Absolute separation histogram")

non_bimodal_likelihoods = []
bimodal_likelihoods = []
chi_non_bimodal = []
chi_bimodal = []

for pickle_file in pickle_files:
    file = open(pickle_file, 'rb')
    Separations, Times, CF_Array_Full, N_sys_total, All_unique_systems, All_unique_systems_L, All_unique_systems_T, Luminosities = pickle.load(file)
    file.close()
    
    CF_median = []
    CF_err = []

    for bit in range(11):
        non_zero_inds = np.where(np.array(CF_Array_Full)[:,bit]>0)[0]
        mean = np.mean(np.array(CF_Array_Full)[:,bit][non_zero_inds])
        std = np.std(np.array(CF_Array_Full)[:,bit][non_zero_inds])
        CF_median.append(mean)
        CF_err.append(std)

    CF_err = np.array(CF_err)
    CF_median = np.array(CF_median)
    
    try:
        non_zero_ind = np.where(CF_per_bin_Tobin>0)[0]
        non_bimodal_likelihood = calc_likelihood(CF_median[non_zero_ind], Abs_median[non_zero_ind], CF_err[non_zero_ind])
        bimodal_likelihood = calc_likelihood(CF_median[non_zero_ind], CF_per_bin_Tobin[non_zero_ind], CF_err[non_zero_ind])
        chi_nb = chi_squared(Abs_median[non_zero_ind], CF_median[non_zero_ind])
        chi_b = chi_squared(CF_per_bin_Tobin[non_zero_ind], CF_median[non_zero_ind])
    except:
        non_zero_ind = np.where(CF_per_bin_Tobin[2:]>0)[0]
        non_bimodal_likelihood = calc_likelihood(CF_median[non_zero_ind], Abs_median[2:][non_zero_ind], CF_err[non_zero_ind])
        bimodal_likelihood = calc_likelihood(CF_median[non_zero_ind], CF_per_bin_Tobin[2:][non_zero_ind], CF_err[non_zero_ind])
        chi_nb = chi_squared(Abs_median[2:][non_zero_ind], CF_median[non_zero_ind])
        chi_b = chi_squared(CF_per_bin_Tobin[2:][non_zero_ind], CF_median[non_zero_ind])
    non_bimodal_likelihoods.append(non_bimodal_likelihood)
    bimodal_likelihoods.append(bimodal_likelihood)
    chi_non_bimodal.append(chi_nb)
    chi_bimodal.append(chi_b)
    print("Updated likelihood for:", pickle_file)

bayes_factor = np.mean(bimodal_likelihood)/np.mean(non_bimodal_likelihood)
print(("Bayes Factor: {0:5.2f} for ".format(bayes_factor)))
