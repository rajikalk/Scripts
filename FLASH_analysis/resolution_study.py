import numpy as np
import pickle
import matplotlib.pyplot as plt
import yt
import sys

r_sink_9 = yt.YTQuantity(4.89593797, 'AU')
L_reference = int(sys.argv[1])
r_scale = (1/(2**(L_reference-9)))*r_sink_9
pickle_files = ['Lref_7.pkl', 'Lref_8.pkl', 'Lref_9.pkl', 'Lref_10.pkl', 'Lref_11.pkl']
R_proto = yt.YTQuantity(2, 'rsun')
r_sink = [4*r_sink_9, 2*r_sink_9, r_sink_9, (1/(2))*r_sink_9, (1/(2**2))*r_sink_9, R_proto.in_units('au')]
#          7              8            9           10                  11 need Lref 18 to resolve protostellar surface
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
line_styles = [':', (0, (3, 5, 1, 5, 1, 5)), '-.', '--', '-']# ['-', '--', '-.', ':']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
plot_it = -1
end_time = 10000
L_scale= np.sqrt(2)/16
L_spec_scale = np.sqrt(2)/2
m_scale = 1/8
scaling_factor = [L_reference-7, L_reference-8, L_reference-9, L_reference-10, L_reference-11]
comp_time = np.nan
comp_h = []
comp_period = []

#============Not scaled semilogy version============================

plot_it = -1
end_time = 10000
L_scale= np.sqrt(2)/16
L_spec_scale = np.sqrt(2)/2
m_scale = 1/8
scaling_factor = [L_reference-7, L_reference-8, L_reference-9, L_reference-10, L_reference-11]
comp_time = np.nan
comp_h = []
comp_period = []

Mass_2500 = []
L_2500 = []
h_2500 = []
T_rot_2500 = []
Mass_err = []
L_err = []
h_err = []
T_rot_err = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    plot_it = plot_it + 1
    file = open(pickle_file, 'rb')
    sink_data, line_counter = pickle.load(file)
    file.close()
    L_ref = int(pickle_file.split('_')[-1].split('.')[0])
    
    form_time = np.inf
    sink_it = -1
    primary_mass = 0
    primary_ind = np.nan
    for sink_id in sink_data.keys():
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        if mass[-1].in_units('msun').value > primary_mass:
            primary_mass = mass[-1].in_units('msun').value
            primary_ind = sink_id
    
    form_time = sink_data[primary_ind]['time'][0]
    if np.isnan(comp_time):
        comp_time = yt.YTQuantity((sink_data[primary_ind]['time'] - form_time)[-1], 's')
    
    for sink_id in sink_data.keys():
        sink_it = sink_it + 1

        time = sink_data[sink_id]['time'] - form_time
        time = yt.YTArray(time, 's')
        end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
        plot_time = time.in_units('yr')[:end_ind+1]
        
        L_spec_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
        L_spec_tot = yt.YTArray(L_spec_tot, 'cm**2/s')
        plot_L_spec_tot = L_spec_tot.in_units('m**2/s')[:end_ind+1]
        
        #Sink_rot = plot_L_spec_tot.in_units('km**2/s') / r_sink[plot_it].in_units('km') #I think this is wrong
        #Sink_circum = 2*np.pi * r_sink[plot_it].in_units('km')
        #Sink_period = Sink_circum/Sink_rot
        
        L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        plot_L_tot = L_tot.in_units('kg*m**2/s')[:end_ind+1]
        
        mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
        plot_mass = mass.in_units('msun')[:end_ind+1]
        
        #===================================================
        #Calculate rotation Period
        r_sink_val = (1/(2**(L_ref-9)))*r_sink_9
        ang_vel = (5/2)*(plot_L_spec_tot.in_units('au**2/day')/(r_sink_val.in_units('au')**2))
        Sink_period = ((4*np.pi)/5) * ((r_sink_val.in_units('au')**2)/(plot_L_spec_tot.in_units('au**2/day')))
        if sink_id == primary_ind:
            ind_2000 = np.argmin(abs(plot_time.value-1500))
            ind_2500 = np.argmin(abs(plot_time.value-2500))
            
            if L_ref != 11:
                
                median_vals = [np.median(plot_L_spec_tot[ind_2000:ind_2500]), np.median(plot_L_tot[ind_2000:ind_2500]), np.median(plot_mass[ind_2000:ind_2500]), np.median(Sink_period[ind_2000:ind_2500])]
                mean_vals = [np.mean(plot_L_spec_tot[ind_2000:ind_2500]), np.mean(plot_L_tot[ind_2000:ind_2500]), np.mean(plot_mass[ind_2000:ind_2500]), np.mean(Sink_period[ind_2000:ind_2500])]
                std_vals = [np.std(plot_L_spec_tot[ind_2000:ind_2500]), np.std(plot_L_tot[ind_2000:ind_2500]), np.std(plot_mass[ind_2000:ind_2500]), np.std(Sink_period[ind_2000:ind_2500])]
                min_vals = [np.min(plot_L_spec_tot[ind_2000:ind_2500]), np.min(plot_L_tot[ind_2000:ind_2500]), np.min(plot_mass[ind_2000:ind_2500]), np.min(Sink_period[ind_2000:ind_2500])]
                max_vals = [np.max(plot_L_spec_tot[ind_2000:ind_2500]), np.max(plot_L_tot[ind_2000:ind_2500]), np.max(plot_mass[ind_2000:ind_2500]), np.max(Sink_period[ind_2000:ind_2500])]
                mid_vals = (np.array(max_vals) + np.array(min_vals))/2
                err_vals = np.array(max_vals) - mid_vals
                
                h_2500.append(yt.YTArray(mid_vals[0], 'm**2/s'))
                L_2500.append(yt.YTArray(mid_vals[1], 'kg*m**2/s'))
                Mass_2500.append(yt.YTArray(mid_vals[2], 'Msun'))
                T_rot_2500.append(yt.YTArray(mid_vals[3], 'day'))
                
                h_err.append(yt.YTArray(err_vals[0], 'm**2/s'))
                L_err.append(yt.YTArray(err_vals[1], 'kg*m**2/s'))
                Mass_err.append(yt.YTArray(err_vals[2], 'Msun'))
                T_rot_err.append(yt.YTArray(err_vals[3], 'day'))
                
                '''
                h_2500.append(mean_vals[0])
                L_2500.append(mean_vals[1])
                Mass_2500.append(mean_vals[2])
                T_rot_2500.append(mean_vals[3])
                
                h_err.append(std_vals[0])
                L_err.append(std_vals[1])
                Mass_err.append(std_vals[2])
                T_rot_err.append(std_vals[3])
                '''
                '''
                h_err.append([median_vals[0] - (mean_vals[0]-std_vals[0]), (mean_vals[0]+std_vals[0])-median_vals[0]])
                L_err.append([median_vals[1] - (mean_vals[1]-std_vals[1]), (mean_vals[1]+std_vals[1])-median_vals[1]])
                Mass_err.append([median_vals[2] - (mean_vals[2]-std_vals[2]), (mean_vals[2]+std_vals[2])-median_vals[2]])
                T_rot_err.append([median_vals[3] - (mean_vals[3]-std_vals[3]), (mean_vals[3]+std_vals[3])-median_vals[3]])
                
                
                h_2500.append(median_vals[0])
                L_2500.append(median_vals[1])
                Mass_2500.append(median_vals[2])
                T_rot_2500.append(median_vals[3])
                
                h_err.append([median_vals[0] - (mean_vals[0]-std_vals[0]), (mean_vals[0]+std_vals[0])-median_vals[0]])
                L_err.append([median_vals[1] - (mean_vals[1]-std_vals[1]), (mean_vals[1]+std_vals[1])-median_vals[1]])
                Mass_err.append([median_vals[2] - (mean_vals[2]-std_vals[2]), (mean_vals[2]+std_vals[2])-median_vals[2]])
                T_rot_err.append([median_vals[3] - (mean_vals[3]-std_vals[3]), (mean_vals[3]+std_vals[3])-median_vals[3]])
                
                Mass_err.append([np.min(plot_mass[ind_2000:ind_2500]), np.max(plot_mass[ind_2000:ind_2500])])
                L_err.append([np.min(plot_L_tot[ind_2000:ind_2500]), np.max(plot_L_tot[ind_2000:ind_2500])])
                h_err.append([np.min(plot_L_spec_tot[ind_2000:ind_2500]), np.max(plot_L_spec_tot[ind_2000:ind_2500])])
                T_rot_err.append([np.min(Sink_period[ind_2000:ind_2500]), np.max(Sink_period[ind_2000:ind_2500])])
                '''
            else:
                h_2500.append(plot_L_spec_tot[-1])
                L_2500.append(plot_L_tot[-1])
                Mass_2500.append(plot_mass[-1])
                T_rot_2500.append(Sink_period[-1])
        
        #================================

    
        if sink_it == 0:
            #if pickle_file == 'Lref_10.pkl':
            #    import pdb
            #    pdb.set_trace()
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, label="$R_{\mathrm{sink}}$="+str(np.round(r_sink[plot_it], decimals=1))+"au", color=colors[plot_it])
            comp_h.append(plot_L_spec_tot[np.argmin(abs(plot_time - comp_time.in_units('yr')))])
            comp_period.append(Sink_period[np.argmin(abs(plot_time - comp_time.in_units('yr')))].in_units('day'))
        else:
            axs.flatten()[0].semilogy(plot_time, plot_mass, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[1].semilogy(plot_time, plot_L_tot, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[2].semilogy(plot_time, plot_L_spec_tot, linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        axs.flatten()[3].semilogy(plot_time, Sink_period.in_units('day'), linestyle=line_styles[plot_it], lw=1, color=colors[plot_it])
        
    
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.3)
axs.flatten()[0].set_ylim([1.e-2, 1.e0])
axs.flatten()[1].set_ylabel("$L$ (kg$\,$m$^2$/s)", labelpad=-0.3)
axs.flatten()[1].set_ylim(bottom=1.e42)
axs.flatten()[2].set_ylabel("$h$ (m$^2$/s)", labelpad=-0.3)
axs.flatten()[2].set_ylim(bottom=1.e14)

axs.flatten()[3].set_ylabel("$P$ (days)", labelpad=-0.3)
axs.flatten()[3].set_xlabel("time (yr)", labelpad=-0.1)
axs.flatten()[3].set_ylim([1.e3, 5.e6])

axs.flatten()[2].set_xlim([0, 10000])
axs.flatten()[0].legend(loc='lower right')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('resolution_study_raw_semilogy.pdf', bbox_inches='tight', pad_inches=0.02)
#plt.savefig('resolution_study_scaled_L'+str(L_reference)+'.pdf', bbox_inches='tight')
        
comp_period.append(yt.YTQuantity(np.nan, 'day'))

plt.clf()
plt.figure(figsize=(single_col_width, 0.3*page_height))
plt.semilogx(r_sink, comp_period)
plt.xlabel('$R_{\mathrm{sink}}$ (au)')
plt.ylabel('$P$ (days)')
plt.savefig("period_vs_r_sink.png", bbox_inches='tight')

#=============Ploting values (2500) vs resolution, with Lref 11
r_sink = [4*r_sink_9, 2*r_sink_9, r_sink_9, (1/(2))*r_sink_9, (1/(2**2))*r_sink_9]

plt.cla()
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

#axs.flatten()[0].loglog(r_sink[:-1], Mass_2500[:-1], lw=1, label='Mass')
axs.flatten()[0].errorbar(r_sink[:-1], np.array(Mass_2500[:-1]), label='Mass', yerr=np.array(Mass_err), fmt='o', marker='.')
axs.flatten()[0].errorbar(r_sink[-1], Mass_2500[-1], yerr=0.25*Mass_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[0].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.01)

#axs.flatten()[1].loglog(r_sink[:-1], L_2500[:-1], lw=1, label='Angular momentum')
axs.flatten()[1].errorbar(r_sink[:-1], np.array(L_2500[:-1]), label='Angular momentum', yerr=np.array(L_err), fmt='o', marker='.')
axs.flatten()[1].errorbar(r_sink[-1], L_2500[-1], yerr=0.7*L_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[1].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[1].set_ylabel("$L$ (kg$\,$m$^2$/s)", labelpad=-0.01)

#axs.flatten()[2].loglog(r_sink[:-1], h_2500[:-1], lw=1, label='Specific Angular momentum')
axs.flatten()[2].errorbar(r_sink[:-1], np.array(h_2500[:-1]), label='Specific Angular momentum', yerr=np.array(h_err), fmt='o', marker='.')
axs.flatten()[2].errorbar(r_sink[-1], h_2500[-1], yerr=0.4*h_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[2].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[2].set_ylabel("$h$ (m$^2$/s)", labelpad=-0.01)

#axs.flatten()[3].loglog(r_sink[:-1], T_rot_2500[:-1], lw=1, label='Rotation period')
axs.flatten()[3].errorbar(r_sink[:-1], np.array(T_rot_2500[:-1]), label='Rotation period', yerr=np.array(T_rot_err), fmt='o', marker='.')
axs.flatten()[3].errorbar(r_sink[-1], T_rot_2500[-1], yerr=0.25*T_rot_2500[-1], uplims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[3].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[3].set_ylabel("$P$ (days)", labelpad=-0.01)

axs.flatten()[3].set_xlabel('$R_\mathrm{sink}$ (au)', labelpad=-0.02)

from scipy.optimize import curve_fit

x_fit = np.linspace(np.log10(R_proto.in_units('au').value), np.log10(r_sink[0].value))
#x_fit = np.linspace(np.log10(r_sink[-1].value), np.log10(r_sink[0].value))
x = 10**x_fit
def power(x, b, c):
    return b*x**c
    
b_err = 0.1
c_err = 0.1

c_guess = (np.max(np.log10(Mass_2500[:-1])) - np.min(np.log10(Mass_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(Mass_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], Mass_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]), sigma=Mass_err[::-1], absolute_sigma=True)
fit_err = np.sqrt(np.diag(pcov))


#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[0].loglog(x, power(x, popt[0], popt[1]), lw=1, color='#ff7f0e')
axs.flatten()[0].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[0].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[0].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='#ff7f0e', alpha=0.15)
M_2rsun = [power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])[0], power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])[0]]
M_eff = [M_2rsun[0]/(Mass_2500[2]+Mass_err[2]).value, M_2rsun[1]/(Mass_2500[2]-Mass_err[2]).value]

M_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
M_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])

eqn_str = '$y=('+str(np.round(popt[0], decimals=2))+'\pm'+str(np.round(fit_err[0], decimals=2))+')x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#eqn_str = "-".join(("".join(eqn_str.split(' + '))).split(' - '))
axs.flatten()[0].text(1.e-2, Mass_2500[1], eqn_str)

b_err = 5
c_err = 5
c_guess = (np.max(np.log10(L_2500[:-1])) - np.min(np.log10(L_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
#b_guess = 10**(np.log10(np.max(L_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
b_guess = np.max(L_2500[:-1])/(np.max(r_sink[:-1])**c_guess)
initial_guess = [b_guess, c_guess]
#popt, pcov = curve_fit(power, r_sink[:-1][::-1], L_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), 0], [10**(np.log10(b_guess)+b_err), c_guess+c_err]), sigma=L_err[::-1], absolute_sigma=True)
popt, pcov = curve_fit(power, r_sink[:-1][::-1], L_2500[:-1][::-1], p0=initial_guess, sigma=np.array(L_err[::-1])/np.array(L_2500[:-1][::-1]), absolute_sigma=False)
fit_err = np.sqrt(np.diag(pcov))

#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[1].loglog(x, power(x, popt[0], popt[1]), lw=1, color='#ff7f0e')
axs.flatten()[1].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[1].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[1].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='#ff7f0e', alpha=0.15)
L_2rsun = [power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])[0], power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])[0]]
L_eff = [L_2rsun[0]/(L_2500[2]+L_err[2]).value, L_2rsun[1]/(L_2500[2]-L_err[2]).value]

L_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
L_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])

#import pdb
#pdb.set_trace()

b_fit_string = '{:.2e}'.format(popt[0]).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')\\times 10^{ '+str(int(b_fit_string[1]))+'})x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#eqn_str = "-".join(("".join(eqn_str.split(' + '))).split(' - '))
axs.flatten()[1].text(1.e-2, L_2500[1], eqn_str)

b_err = 5
c_err = 5
c_guess = (np.max(np.log10(h_2500[:-1])) - np.min(np.log10(h_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(h_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], h_2500[:-1][::-1], p0=initial_guess, sigma=np.array(h_err[::-1])/np.array(h_2500[:-1][::-1]), absolute_sigma=False)
fit_err = np.sqrt(np.diag(pcov))

#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[2].loglog(x, power(x, popt[0], popt[1]), lw=1, color='#ff7f0e')
axs.flatten()[2].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[2].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[2].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='#ff7f0e', alpha=0.15)

h_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
h_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])
#axs.flatten()[2].loglog(x, h_lower, lw=1, c='b')
#axs.flatten()[2].loglog(x, h_upper, lw=1, c='b')
#axs.flatten()[2].fill_between(x, h_lower, h_upper, color='b', alpha=0.25)

h_2rsun = [h_lower[0], h_upper[0]]
h_eff = [h_2rsun[0]/(h_2500[2]+h_err[2]).value, h_2rsun[1]/(h_2500[2]-h_err[2]).value]

b_fit_string = '{:.2e}'.format(popt[0]).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')\\times 10^{ '+str(int(b_fit_string[1]))+'})x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#eqn_str = "-".join(("".join(eqn_str.split(' + '))).split(' - '))

axs.flatten()[2].text(1.e-2, h_2500[1], eqn_str)

c_guess = (np.max(np.log10(T_rot_2500[:-1])) - np.min(np.log10(T_rot_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(T_rot_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
#popt, pcov = curve_fit(power, r_sink[:-1][::-1], T_rot_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]), sigma=T_rot_err[::-1], absolute_sigma=True)
popt, pcov = curve_fit(power, r_sink[:-1][::-1], T_rot_2500[:-1][::-1], p0=initial_guess, sigma=np.array(T_rot_err[::-1])/np.array(T_rot_2500[:-1][::-1]), absolute_sigma=False)
fit_err = np.sqrt(np.diag(pcov))

b_fit_string = '{:.2e}'.format(popt[0]).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')\\times 10^{ '+str(int(b_fit_string[1]))+'})x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[3].loglog(x, power(x, popt[0], popt[1]), lw=1, color='#ff7f0e')
axs.flatten()[3].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[3].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=0.5, color='#ff7f0e')
axs.flatten()[3].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='#ff7f0e', alpha=0.15)

#Sink_period_lower = 1/ang_vel_upper
Sink_period_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
Sink_period_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])

T_2rsun = [Sink_period_lower[0], Sink_period_upper[0]]
T_eff = [T_2rsun[0]/(T_rot_2500[2]+T_rot_err[2]).value, T_2rsun[1]/(T_rot_2500[2]-T_rot_err[2]).value]

#axs.flatten()[3].loglog(x, Sink_period_lower, lw=1, c='b')
#axs.flatten()[3].loglog(x, Sink_period_upper, lw=1, c='b')
#axs.flatten()[3].fill_between(x, Sink_period_lower, Sink_period_upper, color='b', alpha=0.25)

print("Accretion efficiencies of Mass, Angular momentum, Specific Angular momentum, Rotation:", M_eff, L_eff, h_eff, T_eff)
axs.flatten()[3].text(1.e-2, T_rot_2500[1], eqn_str)
axs.flatten()[3].axhline(y=2, c='k', ls='--')

axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('resolution_abs_err.pdf', bbox_inches='tight', pad_inches=0.02)

#=============Ploting values (2500) vs resolution, with Lref 11
r_sink = [4*r_sink_9, 2*r_sink_9, r_sink_9, (1/(2))*r_sink_9, (1/(2**2))*r_sink_9]

plt.cla()
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, 0.7*page_height), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].loglog(r_sink[:-1], Mass_2500[:-1], lw=1, label='Mass')
axs.flatten()[0].errorbar(r_sink[:-1], np.array(Mass_2500[:-1]), label='Mass', yerr=np.array(Mass_err), marker='.')
axs.flatten()[0].errorbar(r_sink[-1], Mass_2500[-1], yerr=0.25*Mass_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[0].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[0].set_ylabel("Mass (M$_\odot$)", labelpad=-0.01)

axs.flatten()[1].loglog(r_sink[:-1], L_2500[:-1], lw=1, label='Angular momentum')
axs.flatten()[1].errorbar(r_sink[:-1], np.array(L_2500[:-1]), label='Angular momentum', yerr=np.array(L_err), marker='.')
axs.flatten()[1].errorbar(r_sink[-1], L_2500[-1], yerr=0.7*L_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[1].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[1].set_ylabel("$L$ (kg$\,$m$^2$/s)", labelpad=-0.01)

axs.flatten()[2].loglog(r_sink[:-1], h_2500[:-1], lw=1, label='Specific Angular momentum')
axs.flatten()[2].errorbar(r_sink[:-1], np.array(h_2500[:-1]), label='Specific Angular momentum', yerr=np.array(h_err), marker='.')
axs.flatten()[2].errorbar(r_sink[-1], h_2500[-1], yerr=0.4*h_2500[-1], lolims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[2].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[2].set_ylabel("$h$ (m$^2$/s)", labelpad=-0.01)

axs.flatten()[3].loglog(r_sink[:-1], T_rot_2500[:-1], lw=1, label='Rotation period')
axs.flatten()[3].errorbar(r_sink[:-1], np.array(T_rot_2500[:-1]), label='Rotation period', yerr=np.array(T_rot_err), marker='.')
axs.flatten()[3].errorbar(r_sink[-1], T_rot_2500[-1], yerr=0.25*T_rot_2500[-1], uplims=[True], marker='.', c='k', alpha=0.25, edgecolor=None)
axs.flatten()[3].axvline(x=R_proto.in_units('au').value, color='k', ls=':', alpha=0.5)
axs.flatten()[3].set_ylabel("$P$ (days)", labelpad=-0.01)

axs.flatten()[3].set_xlabel('$R_\mathrm{sink}$ (au)', labelpad=-0.02)

from scipy.optimize import curve_fit

x_fit = np.linspace(np.log10(R_proto.in_units('au').value), np.log10(r_sink[0].value))
#x_fit = np.linspace(np.log10(r_sink[-1].value), np.log10(r_sink[0].value))
x = 10**x_fit
def power(x, b, c):
    return b*x**c
    
b_err = 0.1
c_err = 0.1

c_guess = (np.max(np.log10(Mass_2500[:-1])) - np.min(np.log10(Mass_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
b_guess = 10**(np.log10(np.max(Mass_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
initial_guess = [b_guess, c_guess]
popt, pcov = curve_fit(power, r_sink[:-1][::-1], Mass_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), c_guess-c_err], [10**(np.log10(b_guess)+b_err), c_guess+c_err]), sigma=Mass_err[::-1], absolute_sigma=True)
fit_err = np.sqrt(np.diag(pcov))


#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[0].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=1, c='b')
axs.flatten()[0].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=1, c='b')
axs.flatten()[0].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='b', alpha=0.25)
M_2rsun = [power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])[0], power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])[0]]
M_eff = [M_2rsun[0]/(Mass_2500[2]+Mass_err[2]).value, M_2rsun[1]/(Mass_2500[2]-Mass_err[2]).value]

M_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
M_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])

eqn_str = '$y=('+str(np.round(popt[0], decimals=2))+'\pm'+str(np.round(fit_err[0], decimals=2))+')x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#eqn_str = "-".join(("".join(eqn_str.split(' + '))).split(' - '))
axs.flatten()[0].text(1.e-2, Mass_2500[1], eqn_str)

b_err = 5
c_err = 5
c_guess = (np.max(np.log10(L_2500[:-1])) - np.min(np.log10(L_2500[:-1])))/(np.max(np.log10(r_sink[:-1])) - np.min(np.log10(r_sink[:-1])))
#b_guess = 10**(np.log10(np.max(L_2500[:-1])) - c_guess*np.log10(np.max(r_sink[:-1])))
b_guess = np.max(L_2500[:-1])/(np.max(r_sink[:-1])**c_guess)
initial_guess = [b_guess, c_guess]
#popt, pcov = curve_fit(power, r_sink[:-1][::-1], L_2500[:-1][::-1], p0=initial_guess, bounds=([10**(np.log10(b_guess)-b_err), 0], [10**(np.log10(b_guess)+b_err), c_guess+c_err]), sigma=L_err[::-1], absolute_sigma=True)
popt, pcov = curve_fit(power, r_sink[:-1][::-1], L_2500[:-1][::-1], p0=initial_guess, sigma=np.array(L_err[::-1])/np.array(L_2500[:-1][::-1]), absolute_sigma=False)
fit_err = np.sqrt(np.diag(pcov))

#axs.flatten()[0].loglog(x, power(x, *popt), lw=1)
axs.flatten()[1].loglog(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), lw=1, c='b')
axs.flatten()[1].loglog(x, power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), lw=1, c='b')
axs.flatten()[1].fill_between(x, power(x, popt[0]+fit_err[0], popt[1]-fit_err[1]), power(x, popt[0]-fit_err[0], popt[1]+fit_err[1]), color='b', alpha=0.25)
L_2rsun = [power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])[0], power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])[0]]
L_eff = [L_2rsun[0]/(L_2500[2]+L_err[2]).value, L_2rsun[1]/(L_2500[2]-L_err[2]).value]

L_lower = power(x, popt[0]-fit_err[0], popt[1]+fit_err[1])
L_upper = power(x, popt[0]+fit_err[0], popt[1]-fit_err[1])

#import pdb
#pdb.set_trace()

b_fit_string = '{:.2e}'.format(popt[0]).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')e'+b_fit_string[1]+')x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'

#eqn_str = "-".join(("".join(eqn_str.split(' + '))).split(' - '))
axs.flatten()[1].text(1.e-2, L_2500[1], eqn_str)

h_lower = L_lower/yt.YTArray(M_upper, 'Msun').in_units('kg')
h_upper = L_upper/yt.YTArray(M_lower, 'Msun').in_units('kg')
axs.flatten()[2].loglog(x, h_lower, lw=1, c='b')
axs.flatten()[2].loglog(x, h_upper, lw=1, c='b')
axs.flatten()[2].fill_between(x, h_lower, h_upper, color='b', alpha=0.25)

h_grad_low = (np.log10(h_lower[-1]) - np.log10(h_lower[0]))/(np.log10(x[-1]) - np.log10(x[0]))
h_grad_high = (np.log10(h_upper[-1]) - np.log10(h_upper[0]))/(np.log10(x[-1]) - np.log10(x[0]))
h_grad = (h_grad_low + h_grad_high)/2
d_h_grad = abs((h_grad_low - h_grad))
h_b_low = np.mean(h_lower/(x**(h_grad_low)))
h_b_high = np.mean(h_upper/(x**(h_grad_high)))
h_b = (h_b_low + h_b_high)/2
d_h_b = abs(h_b-h_b_low)
fit_err = [d_h_b.value, d_h_grad]

h_2rsun = [h_lower[0], h_upper[0]]
h_eff = [h_2rsun[0]/(h_2500[2]+h_err[2]).value, h_2rsun[1]/(h_2500[2]-h_err[2]).value]

b_fit_string = '{:.2e}'.format(h_b.value).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')e'+b_fit_string[1]+')x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'
axs.flatten()[2].text(1.e-2, h_2500[1], eqn_str)

ang_vel_upper = (5/2)*(yt.YTArray(h_upper, 'm**2/s').in_units('au**2/day')/(x**2))
Sink_period_lower = ((4*np.pi)/5) * ((x**2)/(yt.YTArray(h_upper, 'm**2/s').in_units('au**2/day')))

ang_vel_lower = (5/2)*(yt.YTArray(h_lower, 'm**2/s').in_units('au**2/day')/(x**2))
Sink_period_upper = ((4*np.pi)/5) * ((x**2)/(yt.YTArray(h_lower, 'm**2/s').in_units('au**2/day')))

T_grad_low = (np.log10(Sink_period_lower[-1]) - np.log10(Sink_period_lower[0]))/(np.log10(x[-1]) - np.log10(x[0]))
T_grad_high = (np.log10(Sink_period_upper[-1]) - np.log10(Sink_period_upper[0]))/(np.log10(x[-1]) - np.log10(x[0]))
T_grad = (T_grad_low + T_grad_high)/2
d_T_grad = abs((T_grad_low - T_grad_high))
T_b_low = np.mean(Sink_period_lower/(x**(T_grad_low)))
T_b_high = np.mean(Sink_period_upper/(x**(T_grad_high)))
T_b = (T_b_low + T_b_high)/2
d_T_b = abs(T_b-T_b_low)
fit_err = [d_T_b.value, d_T_grad]

T_2rsun = [Sink_period_lower[0], Sink_period_upper[0]]
T_eff = [T_2rsun[0]/(T_rot_2500[2]+T_rot_err[2]).value, T_2rsun[1]/(T_rot_2500[2]-T_rot_err[2]).value]

axs.flatten()[3].loglog(x, Sink_period_lower, lw=1, c='b')
axs.flatten()[3].loglog(x, Sink_period_upper, lw=1, c='b')
axs.flatten()[3].fill_between(x, Sink_period_lower, Sink_period_upper, color='b', alpha=0.25)

print("Accretion efficiencies of Mass, Angular momentum, Specific Angular momentum, Rotation:", M_eff, L_eff, h_eff, T_eff)
b_fit_string = '{:.2e}'.format(T_b.value).split('e')
eqn_str = '$y=(('+b_fit_string[0]+'\pm' + str(np.around(fit_err[0]/(10**(int(b_fit_string[1]))), decimals=2))+')e'+b_fit_string[1]+')x^{'+str(np.round(popt[1], decimals=2))+'\pm'+str(np.round(fit_err[1], decimals=2))+'}$'
axs.flatten()[3].text(1.e-2, T_rot_2500[1], eqn_str)
axs.flatten()[3].axhline(y=2, c='k', ls='--')

axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)

plt.savefig('resolution_extrapolation.pdf', bbox_inches='tight', pad_inches=0.02)
