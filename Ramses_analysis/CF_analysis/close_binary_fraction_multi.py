import numpy as np
import matplotlib.pyplot as plt
import pickle

pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G50.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G100.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G125.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G150.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G200.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G400.pkl']

labels = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, single_col_width))
pit = -1
for pickle_file in pickle_files:
    pit = pit + 1
    file = open(pickle_file, 'rb')
    SFE, Close_Fractions = pickle.load(file)
    file.close()
    
    plt.plot(SFE, Close_Fractions, label=labels[pit])
    
plt.legend(ncol=2)
plt.xlim([0, 0.05])
plt.ylim([0, 1])
plt.xlabel('SFE')
plt.ylabel('Close binary fraction (<100au)')
plt.savefig('close_frac_comp.pdf', bbox_inches='tight', pad_inches=0.02)
    
