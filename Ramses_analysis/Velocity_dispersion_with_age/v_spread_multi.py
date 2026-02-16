import pickle
import numpy as np
import matplotlib.pyplot as plt


pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G50/v_spread.pkl', '/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G100/v_spread.pkl', '/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G125/v_spread.pkl', '/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G150/v_spread.pkl', '/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G200/v_spread.pkl', '/groups/astro/rlk/rlk/Analysis_plots/V_spread_vs_age/G400/v_spread.pkl']

label = ['G50', 'G100', 'G125', 'G150', 'G200', 'G400']

plt.clf()
for pickle_it in range(len(pickle_files)):
    file = open(pickle_files[pickle_it], 'rb')
    V_std_all = pickle.load(file)
    file.close()
    
    plt.plot((np.array(V_std_all).T[0]-np.array(V_std_all).T[0][0]), np.array(V_std_all).T[1], label=label[pickle_it])
    
plt.legend()
plt.xlabel('Time since first sink (yr)')
plt.ylabel('Mean Delta RV')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('sigma_v_vs_time_multi_all.png')

plt.clf()
for pickle_it in range(len(pickle_files)):
    file = open(pickle_files[pickle_it], 'rb')
    V_std_all = pickle.load(file)
    file.close()
    
    plt.plot((np.array(V_std_all).T[0]-np.array(V_std_all).T[0][0]), np.array(V_std_all).T[2], label=label[pickle_it])
    
plt.legend()
plt.xlabel('Time since first sink (yr)')
plt.ylabel('Mean Delta RV')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('sigma_v_vs_time_multi_conv.png')

plt.clf()
for pickle_it in range(len(pickle_files)):
    file = open(pickle_files[pickle_it], 'rb')
    V_std_all = pickle.load(file)
    file.close()
    
    plt.plot((np.array(V_std_all).T[0]-np.array(V_std_all).T[0][0]), np.array(V_std_all).T[3], label=label[pickle_it])
    
plt.legend()
plt.xlabel('Time since first sink (yr)')
plt.ylabel('Mean Delta RV')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('sigma_v_vs_time_multi_intermediate.png')

plt.clf()
for pickle_it in range(len(pickle_files)):
    file = open(pickle_files[pickle_it], 'rb')
    V_std_all = pickle.load(file)
    file.close()
    
    plt.plot((np.array(V_std_all).T[0]-np.array(V_std_all).T[0][0]), np.array(V_std_all).T[4], label=label[pickle_it])
    
plt.legend()
plt.xlabel('Time since first sink (yr)')
plt.ylabel('Mean Delta RV')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('sigma_v_vs_time_multi_hig_mass.png')
    
