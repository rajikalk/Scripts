import numpy as np
import matplotlib.pyplot as plt
import glob
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

def power_law(x, a, k):
    return a * np.power(x, k)
    
def line(x, m, b):
    return m*x + b

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

dirs = ['G50', 'G100', 'G125', 'G150', 'G200', 'G400']
labels=['1500', '3000', '3750', '4500', '6000', '12000']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
line_styles = [':', (0, (3, 1, 1, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), '--', '-']

plt.clf()
plt.figure(figsize=(single_col_width,0.7*single_col_width))
smooth_window = 0.05
for pit in range(len(dirs)):
    file = open(dirs[pit]+'/Power_law_break_30/TPCF.pkl', 'rb')
    SFE, grad, exp_err = pickle.load(file)
    file.close()
    
    SFE = 100*SFE
    grad_smoothed = []
    grad_err_upp = []
    grad_err_low = []
    fitting_error_smoothed = []
    for SFE_it in range(len(SFE)):
        low_SFE = SFE[SFE_it] - smooth_window
        if low_SFE < 0:
            low_SFE = 0
        high_SFE = SFE[SFE_it] + smooth_window
        if high_SFE > 5:
            high_SFE = 5
        low_SFE_it = np.argmin(abs(SFE-low_SFE))
        high_SFE_it = np.argmin(abs(SFE-high_SFE))
        mean_grad = np.mean(grad[low_SFE_it:high_SFE_it])
        median_grad = np.median(grad[low_SFE_it:high_SFE_it])
        std_grad = np.std(grad[low_SFE_it:high_SFE_it])
        err_upper = mean_grad+std_grad
        err_lower = mean_grad-std_grad
        grad_smoothed.append(median_grad)
        grad_err_upp.append(err_upper)
        grad_err_low.append(err_lower)
        smoothed_fit_err = np.mean(exp_err[low_SFE_it:high_SFE_it])
        fitting_error_smoothed.append(smoothed_fit_err)
        
    plt.plot(SFE/100, grad_smoothed, label=labels[pit]+'M$_\odot$', linestyle=line_styles[pit])
    #plt.fill_between(SFE/100, grad_err_low-exp_err, grad_err_upp+exp_err, alpha=0.2)
    plt.fill_between(SFE/100, grad_err_low-np.array(exp_err), grad_err_upp+np.array(exp_err), alpha=0.2)

plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
plt.tick_params(which='both', direction='in')
plt.xlabel('SFE', labelpad=-0.4, fontsize=font_size)
plt.ylabel('TPCF gradient', labelpad=0.0, fontsize=font_size)
plt.xlim([0,0.05])
#plt.ylim([-3.75,-1.25])
plt.legend(loc='best', ncol=2, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
plt.savefig('gradient_comparison.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)
print('plotted gradient comparison')

plt.clf()
plt.figure(figsize=(single_col_width,0.7*single_col_width))
for pit in range(len(dirs)):
    file = open(dirs[pit]+'/Power_law_break_30/SFE_5_TPCF.pkl', 'rb')
    sep_centers, TPCF_frac, TPCF_err, power_law_break_ind, popt1 = pickle.load(file)
    file.close()
    
    if pit == 0:
        ylim = [np.min(TPCF_frac), np.max(TPCF_frac)]
    

    plt.errorbar(10**sep_centers, TPCF_frac, yerr=TPCF_err, fmt = 'o', color=colors[pit], markersize=2)
    plt.loglog(10**sep_centers[:7], 10**line(sep_centers[:7], popt1[0], popt1[1]), linestyle=line_styles[pit], color=colors[pit], label=labels[pit]+'M$_\odot$')
    
plt.tick_params(axis='both', which='major', labelsize=font_size)
plt.tick_params(axis='both', which='minor', labelsize=font_size)
plt.tick_params(which='both', direction='in')
plt.xscale("log")
plt.yscale("log")
plt.ylabel("$1+\\omega(r)$", labelpad=-0.2, fontsize=font_size)
plt.xlabel('Separation (Log$_{10}$(AU))', labelpad=-0.3, fontsize=font_size)
plt.ylim(ylim)
plt.xlim([10, 10**sep_centers[-1]])
plt.legend(loc='upper right', fontsize=font_size, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, borderpad=0.2)
plt.savefig('TPCF.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)
