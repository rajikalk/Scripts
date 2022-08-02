import numpy as np
import matplotlib.pyplot as plt
import yt
import glob
import sys
import collections
import matplotlib as mpl
import pickle
import os
import matplotlib.gridspec as gridspec
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

fit_pickles = ['/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G50/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G100/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G125/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G150/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G200/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G400/chi_squared_fit.pkl']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

plt.clf()
plt.figure(figsize=(single_col_width,0.7*single_col_width))

for fit_pick in range(len(fit_pickles)):
    file = open(fit_pickles[fit_pick], 'rb')
    SFE, reduced_chi_square_tobin, p_values = pickle.load(file)
    file.close()
    
    plt.semilogy(SFE, reduced_chi_square_tobin, label=subplot_titles[fit_pick])
    
plt.legend(loc='best')
plt.xlabel('SFE')
plt.ylabel("Fit (<$\chi^2$>)", labelpad=-0.2)
plt.xlim([0.01, 0.05])
plt.savefig("reduced_chi_squared_multi.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02)
