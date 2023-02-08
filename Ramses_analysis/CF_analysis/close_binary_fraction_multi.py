import numpy as np
import matplotlib.pyplot as plt
import pickle
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

#pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G50/G50.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G100/G100.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G125/G125.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G150/G150.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G200/G200.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G400/G400.pkl']

pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G50/G50.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G100/G100.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G125/G125.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G150/G150.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G200/G200.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Close_binary_fraction/G400/G400.pkl']

labels = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
line_styles = [':', (0, (3, 1, 1, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), '--', '-']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.7*single_col_width))
pit = -1
smooth_window = 0.05
for pickle_file in pickle_files:
    pit = pit + 1
    file = open(pickle_file, 'rb')
    SFE, Close_Fractions = pickle.load(file)
    file.close()
    
    SFE = 100*SFE
    frac_smoothed = []
    frac_err_upp = []
    frac_err_low = []
    for SFE_it in range(len(SFE)):
        low_SFE = SFE[SFE_it] - smooth_window
        if low_SFE < 0:
            low_SFE = 0
        high_SFE = SFE[SFE_it] + smooth_window
        if high_SFE > 5:
            high_SFE = 5
        low_SFE_it = np.argmin(abs(SFE-low_SFE))
        high_SFE_it = np.argmin(abs(SFE-high_SFE))
        mean_frac = np.mean(Close_Fractions[low_SFE_it:high_SFE_it])
        median_frac = np.median(Close_Fractions[low_SFE_it:high_SFE_it])
        std_frac = np.std(Close_Fractions[low_SFE_it:high_SFE_it])
        err_upper = mean_frac+std_frac
        err_lower = mean_frac-std_frac
        frac_smoothed.append(median_frac)
        frac_err_upp.append(err_upper)
        frac_err_low.append(err_lower)
    
    plt.plot(SFE, Close_Fractions, label=labels[pit], color=colors[pit], linestyle=line_styles[pit])
    
plt.legend(ncol=2, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
plt.xlim([0, 0.05])
plt.ylim([0, 1])
plt.xlabel('SFE')
plt.ylabel('Close binary fraction (<100au)')
plt.savefig('close_frac_comp.pdf', bbox_inches='tight', pad_inches=0.02)
    
