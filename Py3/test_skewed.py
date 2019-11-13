import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math as math
import scipy.special as sp
import pickle
import matplotlib.gridspec as gridspec

def func(x, sigma, mu, alpha, c, amp):
    #normal distribution
    normpdf = (1/(sigma*np.sqrt(2*math.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigma,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigma))/(np.sqrt(2)))))
    return 2*amp*normpdf*normcdf + c

files = ["Mach_0.1/multiple_folds_over_5_orbits.pkl","Mach_0.2/multiple_folds_over_5_orbits.pkl"]
max_accretion = []
base_accretion = []
strength = []
beta = []
y_fits = []
plot_e = []
for file in files:
    file_open = open(file, 'rb')
    multiple_folds, phase_centers, mean_eccentricity, std_eccentricity, accretion_err, n_lines, multiple_folds_normalised = pickle.load(file_open)
    file_open.close()

    plot_e.append(mean_eccentricity)
    x_data = phase_centers[23:-15]
    x = np.linspace(np.min(x_data),np.max(x_data),100)
    max_accretion.append([])
    base_accretion.append([])
    beta.append([])
    strength.append([])
    y_fits.append([])
    file_name = file.split('/')[0] +'/'
    for orbit in range(len(multiple_folds_normalised)):
        '''
        if mean_eccentricity[orbit] == 0.27:
            import pdb
            pdb.set_trace()
        '''
        y_data = multiple_folds_normalised[orbit][23:-15]
        
        plt.clf()
        plt.plot(x_data,y_data,ls='steps-mid')
        results = []
        for tries in range(50):
            sigma = np.random.random()*2*0.15
            amp = np.random.random()*2*np.max(y_data)
            p = np.array([sigma, x_data[np.argmax(y_data)], -5,np.min(y_data),amp])
            try:
                popt, pcov = curve_fit(func, x_data, y_data, p)
            except:
                pass
            err = np.sum(np.abs(func(x_data, *popt) - y_data))
            results.append((err, popt))
            if err < 0.1:
                break
        err, popt = min(results, key=lambda x:x[0])
        if mean_eccentricity[orbit] == 0.27:
            popt = np.array([0.35, x_data[np.argmax(y_data)]+0.15, -5,np.median(y_data)-0.5,np.max(y_data)*0.2])
        y_fit= func(x, *popt)
        sigmag, mu, alpha, base, amp = popt
        max = np.max(y_fit)
        max_accretion[-1].append(max)
        base_accretion[-1].append(np.min(y_fit))
        beta[-1].append(max/np.min(y_fit))
        strength[-1].append(sigmag)
        plt.plot(x,y_fit)
        plt.ylim([0,6])
        y_fits[-1].append(y_fit)
        print('---------------------------------------------')
        print('eccentricity   = '+str(mean_eccentricity[orbit]))
        print('amplitude      = '+str(amp))
        print('maximum_value  = '+str(np.max(y_fit)))
        print('base_accretion = '+str(np.min(y_fit)))
        print('strength       = '+str(sigmag))
        plt.savefig(file_name+'fittted_eccentricity_'+str(mean_eccentricity[orbit])+'.pdf')
        
#Make normalised fits plot
plt.clf()
fig = plt.figure()
fig.set_size_inches(4.0, 6.0)
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.0)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0], sharex=ax1, sharey=ax1)

n_lines = len(y_fits[0])
c_index = np.linspace(0.0, 0.95, n_lines)

e_int = 0
for fit in y_fits[0]:
    ax1.plot(x, fit, color=plt.cm.magma(c_index[e_int]), label='e='+str(plot_e[0][e_int]))
    e_int = e_int + 1

ax1.legend(loc='center left', bbox_to_anchor=(0.985, 0.5))
ax1.set_ylabel("Normalised Accretion")
xticklabels = ax1.get_xticklabels()
plt.setp(xticklabels, visible=False)
ax1.tick_params(axis='x', which='major', direction="in")

e_int = 0
for fit in y_fits[1]:
    ax2.plot(x, fit, color=plt.cm.magma(c_index[e_int]), label='e='+str(plot_e[0][e_int]))
    e_int = e_int + 1

ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax2.set_xlabel("Orbital Phase ($\phi$)")
ax2.set_ylabel("Normalised Accretion")
ax2.text(0.1, ax1.get_ylim()[1]*0.9, 'T2', va="center", ha="left", color='k', fontsize=args.text_font)
ax1.text(0.1, ax1.get_ylim()[1]*0.9, 'T1', va="center", ha="left", color='k', fontsize=args.text_font)
ax2.tick_params(axis='x', which='major', direction="in")
yticklabels = ax2.get_yticklabels()
plt.setp(yticklabels[-1], visible=False)

plt.savefig('normalised_fits.eps', bbox_inches='tight', pad_inches = 0.02)
plt.savefig('normalised_fits.pdf', bbox_inches='tight', pad_inches = 0.02)

#make beta plot
plt.clf()
plt.scatter(mean_eccentricity, beta[0], label='T1', marker='o')
plt.scatter(mean_eccentricity, beta[1], label='T2', marker='^')
plt.xlabel('eccentricity')
plt.ylabel('$\\beta$')
plt.legend(loc='best')
plt.savefig('beta.pdf')
