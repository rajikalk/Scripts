import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from scipy.optimize import curve_fit

def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Gaussian_cdf(x,scale,mean,sigma):
    return scale*stats.norm.cdf(x, mean, sigma)

def Skewed_Gaussian(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.pdf(x, skew, loc=mean, scale=sigma)
    
def Skewed_Gaussian_cdf(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.cdf(x, skew, loc=mean, scale=sigma)

grad_pickles = ['/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G50/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G100/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G125/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G150/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G200/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G400/grad_pickle.pkl']

Sim_ids = ['G50', 'G100', 'G125', 'G150', 'G200', 'G400']
#Defining gradient bins and getting tick labels
grad_bins = np.concatenate((-1*np.logspace(2,-5,15), np.array([0, 1.e10])))#np.concatenate((-1*np.logspace(2,-6,28)[1:-2], np.array([0, 1.e10])))#np.concatenate((-1*np.logspace(3,-6,19)[1:-2], np.array([0, 1.e10]))) #np.concatenate((-1*np.logspace(5,-3,9), np.array([0, 1.e10])))
#inpiral only:
grad_bin_centers = (grad_bins[1:] + grad_bins[:-1])/2

bin_widths = grad_bins[1:] - grad_bins[:-1]
ticklabels = []
for bin_val in grad_bins:
    if int(np.sign(bin_val)) < 0:
        tick_str = str(int(np.sign(bin_val)))+'0$^{'+str(int(np.log10(abs(bin_val))))+'}$'
    elif int(np.sign(bin_val)) > 0:
        tick_str = '>0'
    else:
        tick_str = '0'
    ticklabels.append(r'{}'.format(tick_str))
ticklabels.append("")

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig1, axs1 = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
fig2, axs2 = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
fig3, axs3 = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
fig4, axs4 = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
fig_list = [fig1, fig2, fig3, fig4]
axs_list = [axs1, axs2, axs3, axs4]
iter_range = range(0, len(grad_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for grad_it in range(len(grad_pickles)):
    file = open(grad_pickles[grad_it], 'rb')
    Initial_gradients, Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000, Grad_1e4, Grad_1e3 = pickle.load(file)
    file.close()

    #Plotting initial gradients
    #grad_bins = np.concatenate((-1*np.logspace(4,-3,15), np.array([0, 1.e10])))
    grad_hist_core, grad_bins = np.histogram(Initial_gradients[0], bins=grad_bins)
    grad_hist_core_delayed, grad_bins = np.histogram(Initial_gradients[1], bins=grad_bins)
    grad_hist_capt, grad_bins = np.histogram(Initial_gradients[2], bins=grad_bins)
    grad_hist_misc, grad_bins = np.histogram(Initial_gradients[3], bins=grad_bins)

    x_range = np.arange(len(grad_hist_core)+1)

    grad_hist_core_norm = grad_hist_core/np.sum(grad_hist_core)
    grad_hist_core_delayed_norm = grad_hist_core_delayed/np.sum(grad_hist_core_delayed)
    grad_hist_capt_norm = grad_hist_capt/np.sum(grad_hist_capt)
    grad_hist_misc_norm = grad_hist_misc/np.sum(grad_hist_misc)
    
    grad_hist_core_norm = np.concatenate((grad_hist_core_norm, np.array([grad_hist_core_norm[-1]])))
    grad_hist_core_delayed_norm = np.concatenate((grad_hist_core_delayed_norm, np.array([grad_hist_core_delayed_norm[-1]])))
    grad_hist_capt_norm = np.concatenate((grad_hist_capt_norm, np.array([grad_hist_capt_norm[-1]])))
    grad_hist_misc_norm = np.concatenate((grad_hist_misc_norm, np.array([grad_hist_misc_norm[-1]])))

    axs_list[0][grad_it].step(x_range, grad_hist_core_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
    axs_list[0][grad_it].step(x_range, grad_hist_core_delayed_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
    axs_list[0][grad_it].step(x_range, grad_hist_capt_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
    axs_list[0][grad_it].step(x_range, grad_hist_misc_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
    
    axs_list[0][grad_it].set_ylabel('')
    
    axs_list[0][grad_it].tick_params(axis='both', which='major', labelsize=font_size)
    axs_list[0][grad_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs_list[0][grad_it].tick_params(axis='x', direction='in')
    axs_list[0][grad_it].tick_params(axis='y', direction='in')
    '''
    if grad_it == 0:
        axs1[grad_it].legend(loc='best')
    '''
    '''
    if grad_it == (len(grad_pickles) - 1):
        axs_list[0][grad_it].set_xlim([x_range[0], x_range[-1]])
        axs_list[0][grad_it].set_xticks(x_range[::2])
        axs_list[0][grad_it].set_xticklabels(ticklabels[::2][:-1])
        axs_list[0][grad_it].set_xlabel('Inspiral rate (au/yr)')
        axs_list[0][grad_it].set_ylim(bottom=0)
    '''
    fig_list[0].savefig('Initial_grad_hist.png', bbox_inches='tight', pad_inches=0.02)

    mean_grads = [Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000]

    #Plotting mean gradients
    #calculate means
    time_means = [1000, 10000, 100000]
    time_means_counter = 0
    for Initial_mean_grad in mean_grads:
        plt.subplots_adjust(hspace=-0.1)
        core_mean = []
        for grit in range(len(Initial_mean_grad[0])):
            core_mean.append(np.mean(Initial_mean_grad[0][grit]))
        grad_hist_core_mean, grad_bins = np.histogram(core_mean, bins=grad_bins)
        grad_hist_core_mean_std = np.sqrt(grad_hist_core_mean)
        grad_hist_core_mean_rel_err = grad_hist_core_mean_std/grad_hist_core_mean
        grad_hist_core_mean_rel_err = np.nan_to_num(grad_hist_core_mean_rel_err)
        #grad_hist_core_mean_err = 1/np.sqrt(len(core_mean))
        #grad_hist_core_mean_rel_err = grad_hist_core_mean_err/grad_hist_core_mean
        #grad_hist_core_mean_rel_err[np.isinf(grad_hist_core_mean_rel_err)] = 0

        core_delayed_mean = []
        for grit in range(len(Initial_mean_grad[1])):
            core_delayed_mean.append(np.mean(Initial_mean_grad[1][grit]))
        grad_hist_core_delayed_mean, grad_bins = np.histogram(core_delayed_mean, bins=grad_bins)
        grad_hist_core_delayed_mean_std = np.sqrt(grad_hist_core_delayed_mean)
        grad_hist_core_delayed_mean_rel_err = grad_hist_core_delayed_mean_std/grad_hist_core_delayed_mean
        grad_hist_core_delayed_mean_rel_err = np.nan_to_num(grad_hist_core_delayed_mean_rel_err)
        #grad_hist_core_delayed_mean_err = 1/np.sqrt(len(core_delayed_mean))
        #grad_hist_core_delayed_mean_rel_err = grad_hist_core_delayed_mean_err/grad_hist_core_delayed_mean
        #grad_hist_core_delayed_mean_rel_err[np.isinf(grad_hist_core_delayed_mean_rel_err)] = 0

        capt_mean = []
        for grit in range(len(Initial_mean_grad[2])):
            capt_mean.append(np.mean(Initial_mean_grad[2][grit]))
        grad_hist_capt_mean, grad_bins = np.histogram(capt_mean, bins=grad_bins)
        grad_hist_capt_mean_std = np.sqrt(grad_hist_capt_mean)
        grad_hist_capt_mean_rel_err = grad_hist_capt_mean_std/grad_hist_capt_mean
        grad_hist_capt_mean_rel_err = np.nan_to_num(grad_hist_capt_mean_rel_err)

        misc_mean = []
        for grit in range(len(Initial_mean_grad[3])):
            misc_mean.append(np.mean(Initial_mean_grad[3][grit]))
        grad_hist_misc_mean, grad_bins = np.histogram(misc_mean, bins=grad_bins)
        grad_hist_misc_mean_std = np.sqrt(grad_hist_misc_mean)
        grad_hist_misc_mean_rel_err = grad_hist_misc_mean_std/grad_hist_misc_mean
        grad_hist_misc_mean_rel_err = np.nan_to_num(grad_hist_misc_mean_rel_err)

        grad_hist_core_mean_norm = grad_hist_core_mean/np.sum(grad_hist_core_mean)
        grad_hist_core_delayed_mean_norm = grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean)
        grad_hist_capt_mean_norm = grad_hist_capt_mean/np.sum(grad_hist_capt_mean)
        grad_hist_misc_mean_norm = grad_hist_misc_mean/np.sum(grad_hist_misc_mean)

        grad_hist_core_mean_norm = np.concatenate((grad_hist_core_mean_norm, np.array([grad_hist_core_mean_norm[-1]])))
        grad_hist_core_delayed_mean_norm = np.concatenate((grad_hist_core_delayed_mean_norm, np.array([grad_hist_core_delayed_mean_norm[-1]])))
        grad_hist_capt_mean_norm = np.concatenate((grad_hist_capt_mean_norm, np.array([grad_hist_capt_mean_norm[-1]])))
        grad_hist_misc_mean_norm = np.concatenate((grad_hist_misc_mean_norm, np.array([grad_hist_misc_mean_norm[-1]])))

        grad_hist_core_mean_rel_err = np.concatenate((grad_hist_core_mean_rel_err, np.array([grad_hist_core_mean_rel_err[-1]])))
        grad_hist_core_delayed_mean_rel_err = np.concatenate((grad_hist_core_delayed_mean_rel_err, np.array([grad_hist_core_delayed_mean_rel_err[-1]])))
        grad_hist_capt_mean_rel_err = np.concatenate((grad_hist_capt_mean_rel_err, np.array([grad_hist_capt_mean_rel_err[-1]])))
        grad_hist_misc_mean_rel_err = np.concatenate((grad_hist_misc_mean_rel_err, np.array([grad_hist_misc_mean_rel_err[-1]])))

        axs_list[time_means_counter+1][grad_it].step(x_range, grad_hist_core_mean_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
        if time_means[time_means_counter] == 10000:
            scale_guess = np.max(grad_hist_core_mean_norm[:-1])
            import pdb
            pdb.set_trace()
            mean_guess = np.nanmean(core_mean)
            std_guess = np.nanstd(core_mean)
            popt, pcov = curve_fit(Gaussian, x_range[:-1], (grad_hist_core_mean_norm[:-1]), [scale_guess, mean_guess, std_guess])
            x_fit = np.linspace(0, len(grad_hist_core))
            fit = Gaussian(x_fit, *popt)
            axs_list[time_means_counter+1][grad_it].plot(x_fit, fit, color='b')
        #axs_list[time_means_counter+1][grad_it].errorbar(x_range+0.5, grad_hist_core_mean_norm, yerr=(grad_hist_core_mean_rel_err*grad_hist_core_mean_norm), fmt='none', linewidth=2, color='b', alpha=0.5)
        axs_list[time_means_counter+1][grad_it].step(x_range, grad_hist_core_delayed_mean_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
        if time_means[time_means_counter] == 10000:
            scale_guess = np.max(grad_hist_core_delayed_mean_norm[:-1])
            mean_guess = np.nanmean(core_delayed_mean)
            std_guess = np.nanstd(core_delayed_mean)
            popt, pcov = curve_fit(Gaussian, x_range[:-1], (grad_hist_core_delayed_mean_norm[:-1]), [scale_guess, mean_guess, std_guess])
            x_fit = np.linspace(0, len(grad_hist_core))
            fit = Gaussian(x_fit, *popt)
            axs_list[time_means_counter+1][grad_it].plot(x_fit, fit, color='purple')
        
        #axs_list[time_means_counter+1][grad_it].errorbar(x_range+0.5, grad_hist_core_delayed_mean_norm, yerr=(grad_hist_core_delayed_mean_rel_err*grad_hist_core_delayed_mean_norm), fmt='none', linewidth=2, color='purple', alpha=0.5)
        axs_list[time_means_counter+1][grad_it].step(x_range, grad_hist_capt_mean_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
        if time_means[time_means_counter] == 10000:
            scale_guess = np.max(grad_hist_capt_mean_norm[:-1])
            mean_guess = np.nanmean(capt_mean)
            std_guess = np.nanstd(capt_mean)
            try:
                popt, pcov = curve_fit(Gaussian, x_range[:-1], (grad_hist_capt_mean_norm[:-1]), [scale_guess, mean_guess, std_guess])
                x_fit = np.linspace(0, len(grad_hist_core))
                fit = Gaussian(x_fit, *popt)
                axs_list[time_means_counter+1][grad_it].plot(x_fit, fit, color='red')
            except:
                print("No fit found for Dynamical Capture")
        
        #axs_list[time_means_counter+1][grad_it].errorbar(x_range+0.5, grad_hist_capt_mean_norm, yerr=(grad_hist_capt_mean_rel_err*grad_hist_capt_mean_norm), fmt='none', linewidth=2, color='red', alpha=0.5)
        axs_list[time_means_counter+1][grad_it].step(x_range, grad_hist_misc_mean_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
        if time_means[time_means_counter] == 10000:
            scale_guess = np.max(grad_hist_misc_mean_norm[:-1])
            mean_guess = np.nanmean(misc_mean)
            std_guess = np.nanstd(misc_mean)
            try:
                popt, pcov = curve_fit(Gaussian, x_range[:-1], (grad_hist_misc_mean_norm[:-1]), [scale_guess, mean_guess, std_guess])
                x_fit = np.linspace(0, len(grad_hist_core))
                fit = Gaussian(x_fit, *popt)
                axs_list[time_means_counter+1][grad_it].plot(x_fit, fit, color='orange')
            except:
                print("No fit found for Misc")
        #axs_list[time_means_counter+1][grad_it].errorbar(x_range+0.5, grad_hist_misc_mean_norm, yerr=(grad_hist_misc_mean_rel_err*grad_hist_misc_mean_norm), fmt='none', linewidth=2, color='orange', alpha=0.5)
        axs_list[time_means_counter+1][grad_it].set_xlim([x_range[0], x_range[-1]])
        
        axs_list[time_means_counter+1][grad_it].tick_params(axis='both', which='major', labelsize=font_size)
        axs_list[time_means_counter+1][grad_it].tick_params(axis='both', which='minor', labelsize=font_size)
        axs_list[time_means_counter+1][grad_it].tick_params(axis='x', direction='in')
        axs_list[time_means_counter+1][grad_it].tick_params(axis='y', direction='in')
        '''
        ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_mean/np.sum(grad_hist_core_mean), label="Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='b')
        ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean), label="Delayed Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='purple')
        ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_capt_mean/np.sum(grad_hist_capt_mean), label="Dynamical Capture", width=1, color='None', linewidth=2, edgecolor='red')
        ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_misc_mean/np.sum(grad_hist_misc_mean), label="Other", width=1, color='None', linewidth=2, edgecolor='orange')
        '''
        axs_list[time_means_counter+1][grad_it].set_ylabel('#')
        '''
        if grad_it == 0:
            axs1[time_means_counter+1].legend(loc='best')
        '''
        '''
        if grad_it == (len(grad_pickles) - 1):
            axs_list[time_means_counter+1][grad_it].set_xlim([x_range[0], x_range[-1]])
            axs_list[time_means_counter+1][grad_it].set_xticks(x_range[::2])
            axs_list[time_means_counter+1][grad_it].set_xticklabels(ticklabels[::2][:-1])
            axs_list[time_means_counter+1][grad_it].set_xlabel('Semi-major axis inspiral rate (au/yr)')
            axs_list[time_means_counter+1][grad_it].set_ylim(bottom=0)
        '''
        axs_list[time_means_counter+1][grad_it].set_ylim(bottom=0)
        
        ticklabels = []
        for bin_val in grad_bins:
            if int(np.sign(bin_val)) < 0:
                tick_str = str(int(np.sign(bin_val)))+'0$^{'+str(int(np.log10(abs(bin_val))))+'}$'
            elif int(np.sign(bin_val)) > 0:
                tick_str = '>0'
            else:
                tick_str = '0'
            ticklabels.append(r'{}'.format(tick_str))
        ticklabels.append("")
        
        axs_list[time_means_counter+1][grad_it].set_xticklabels(ticklabels[::2])
        if grad_it == 5:
            axs_list[time_means_counter+1][grad_it].set_xlabel('Inspiral rate (au/yr)')
        axs_list[time_means_counter+1][grad_it].set_ylabel('#')
        
        fig_list[time_means_counter+1].savefig('Initial_mean_grad_'+str(time_means[time_means_counter])+'.png', bbox_inches='tight', pad_inches=0.02)
        time_means_counter = time_means_counter + 1

"""

def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Skewed_Gaussian(x, scale, mean, sigma, skew, kurt):
    return scale*stats.skewnorm.pdf(x,mean, sigma, skew, kurt)
    
#Fitting
import scipy.stats as stats
from scipy.optimize import curve_fit
x_log = np.log10(-1*grad_bins[:-2])-0.5
x_log_fit = np.linspace(x_log[0], x_log[-1], 100)

y_log_core = grad_hist_core_mean_norm[:-2]#[::-1]
y_log_core_err = (grad_hist_core_mean_rel_err*grad_hist_core_mean_norm)[:-2]#[::-1]
y_log_core_delayed = grad_hist_core_delayed_mean_norm[:-2]#[::-1]
y_log_core_delayed_err = (grad_hist_core_delayed_mean_rel_err*grad_hist_core_delayed_mean_norm)[:-2]#[::-1]
y_log_capt = grad_hist_capt_mean_norm[:-2]#[::-1]
y_log_capt_err = (grad_hist_capt_mean_rel_err*grad_hist_capt_mean_norm)[:-2]#[::-1]
y_log_misc = grad_hist_misc_mean_norm[:-2]#[::-1]
y_log_misc_err = (grad_hist_misc_mean_rel_err*grad_hist_misc_mean_norm)[:-2]#[::-1]

mean_guess = x_log[np.argmax(y_log_core)]#-0.5
mode_guess = x_log[np.argmax(y_log_core)]

sig_guess = np.std(np.log10(-1*np.array(core_mean)[np.array(core_mean)<0]))
scale_guess = np.sum(grad_hist_core_mean[:-2])/np.sum(grad_hist_core_mean)
popt_core, pcov_core = curve_fit(Gaussian, x_log, y_log_core, [scale_guess, mean_guess, sig_guess], sigma=y_log_core_err)#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
fit_core = Gaussian(x_log_fit, *popt_core)
print('Core frag mean inspiral:', popt_core[1], '\pm', popt_core[2])

usable_inds = np.where((np.isnan(core_delayed_mean)==False)&(np.array(core_delayed_mean)<0))[0]
core_delayed_power = np.log10(-1*np.array(core_delayed_mean)[usable_inds])
mean_guess = np.mean(core_delayed_power)#-0.5
median_guess = np.median(core_delayed_power)
mode_guess = x_log[np.argmax(y_log_core_delayed)]
sig_guess = np.std(core_delayed_power)
skew_guess = (3*(mean_guess - median_guess))/sig_guess
scale_guess = np.sum(grad_hist_core_delayed_mean[:-2])/np.sum(grad_hist_core_delayed_mean)
popt_core_delayed, pcov_core_delayed = curve_fit(Gaussian, x_log, y_log_core_delayed, [scale_guess, mean_guess, sig_guess], sigma=y_log_core_delayed_err)
fit_core_delayed = Gaussian(x_log_fit, *popt_core_delayed)
print('Delayed core frag mean inspiral:', popt_core_delayed[1], '\pm', popt_core_delayed[2])

plt.clf()
plt.errorbar(x_log, y_log_core_delayed, yerr=y_log_core_delayed_err, marker='o', color='purple', linestyle="")
#y_skew = scale_guess*skewnorm.pdf(x_log_fit, skew_guess, loc=mode, scale=sig_guess)
y_skew = skewnorm.pdf(x_log_fit, -4, loc=2, scale=2)
plt.plot(x_log_fit, y_skew)
plt.xlim([x_log[0], x_log[-1]])
plt.ylim(bottom=0)
plt.savefig('skew_test.png')

mean_guess = x_log[np.argmax(y_log_capt)]#-0.5
sig_guess = np.std(np.log10(-1*np.array(capt_mean)[np.array(capt_mean)<0]))
scale_guess = np.sum(grad_hist_capt_mean[:-2])/np.sum(grad_hist_capt_mean)
popt_capt, pcov_core = curve_fit(Gaussian, x_log, y_log_capt, [scale_guess, mean_guess, sig_guess], sigma=y_log_capt_err)
fit_capt = Gaussian(x_log_fit, *popt_capt)
print('Core frag mean inspiral:', popt_capt[1], '\pm', popt_capt[2])

mean_guess = x_log[np.argmax(y_log_misc)]#-0.5
sig_guess = np.std(np.log10(-1*np.array(misc_mean)[np.array(misc_mean)<0]))
scale_guess = np.sum(grad_hist_misc_mean[:-2])/np.sum(grad_hist_misc_mean)
popt_misc, pcov_misc = curve_fit(Gaussian, x_log, y_log_misc, [scale_guess, mean_guess, sig_guess], sigma=y_log_misc_err)
fit_misc = Gaussian(x_log_fit, *popt_misc)
print('Core frag mean inspiral:', popt_misc[1], '\pm', popt_misc[2])


plt.clf()
plt.errorbar(x_log, y_log_core, yerr=y_log_core_err, marker='o', color='b', linestyle="")
plt.plot(x_log_fit, fit_core, color='b')

plt.errorbar(x_log, y_log_core_delayed, yerr=y_log_core_delayed_err, marker='o', color='purple', linestyle="")
plt.plot(x_log_fit, fit_core_delayed, color='purple')

plt.errorbar(x_log, y_log_capt, yerr=y_log_capt_err, marker='o', color='red', linestyle="")
plt.plot(x_log_fit, fit_capt, color='red')

plt.errorbar(x_log, y_log_misc, yerr=y_log_misc_err, marker='o', color='orange', linestyle="")
plt.plot(x_log_fit, fit_misc, color='orange')
#plt.xscale("log", nonpositive='clip')
#plt.semilogx(x_log, y_log)
plt.xlim([x_log[0], x_log[-1]])
plt.ylim(bottom=0)
plt.savefig('log_normal_fit.png')
"""
