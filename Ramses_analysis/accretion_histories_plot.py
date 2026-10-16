import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib
import matplotlib.pyplot as plt
import pickle

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}" "\sisetup{detect-all}" r"\usepackage{helvet}" r"\usepackage{sansmath}" "\sansmath"

particle_pickle = "/g/data/ek9/rlk100/RAMSES/Global/stars_red_512.pkl"
file = open(particle_pickle, 'rb')
particle_data = pickle.load(file)
file.close()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

#Plot sink histories for classification
plot_individual = False
if plot_individual == False:
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width*0.5))#, sharey=True)

t_0 = particle_data['time'].T[0][0]
for sink_id in range(len(particle_data['dm'][0])):
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-t_0) * units['time_unit'].in_units('Myr')
    if plot_individual == True:
        plt.clf()
        plt.semilogy(time_arr, acc_val)
        plt.xlabel("time")
        plt.ylabel("Accretion rate")
        plt.ylim([1.e-9, 5.e-5])
        plt.title("Sink ID = "+str(sink_id))
        plt.savefig("Accretion_history_"+str(sink_id)+".png")
        print("saved accretion history for sink", sink_id)
    else:
        axs.semilogy(time_arr, acc_val, linewidth=0.25, alpha=0.05)
        if sink_id == 0:
            axs.set_xlim([time_arr[0], time_arr[-1]])
if plot_individual == False:
    axs.set_xlabel("Simulation time (Myr)")
    axs.set_ylabel("Accretion rate ($\mathrm{M}_\odot/\mathrm{yr}$)")
    axs.set_ylim([1.e-9, 5.e-5])
    plt.savefig("All_accretion_histories.pdf", bbox_inches='tight', pad_inches=0.02)
    print("saved all accretion histories")

IDs_continuous_clean = [2, 8, 11, 89, 90, 126, 129, 162, 165, 208,  222, 219, 218]
IDs_continuous_messy = [84, 109, 127, 139, 140, 141, 153, 157, 159, 164, 166, 170, 172, 184, 185, 214, 221]

IDs_suddenly_shut_off_clean = [0, 1, 59, 64, 96, 115, 117, 167, 196, 217]
IDs_suddenly_shut_off_messy = [17, 20, 26, 49, 110, 147, 151, 154, 173, 175, 181]

IDs_steady_decline_clean = [5, 12, 40, 42, 57, 65, 95, 104, 119, 112, 130, 136, 142, 145, 158, 163, 169, 171, 194]
IDs_steady_decline_messy = [18, 19, 25, 28, 35, 38, 46, 51, 66, 69, 97, 116, 146, 188, 192, 200]

IDs_late_accretion_clean = [6, 14, 27, 79, 82, 100, 124, 125, 135, 138, 149, 179, 180, 220]
IDs_late_accretion_messy = [13, 24, 30, 62, 112, 144, 156, 160, 168, 178, 203, 216]

IDs_Other = [3, 9, 10, 22, 31, 33, 37, 44, 48, 50, 55, 92, 93, 94, 98, 111, 114, 118, 123, 128, 143, 177, 182, 187, 193, 195, 197, 207, 201, 210, 223]

ln_wd = 0.5
alpha = 1

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, single_col_width), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0.0)

for sink_id in IDs_continuous_clean:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[0].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)

axs.flatten()[0].set_ylim([1.e-9, 5.e-5])
axs.flatten()[0].set_xlim(left=0.0)

'''
for sink_id in IDs_continuous_messy:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[0].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)
'''
for sink_id in IDs_suddenly_shut_off_clean:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[1].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)
'''
for sink_id in IDs_suddenly_shut_off_messy:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[1].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)
'''
for sink_id in IDs_steady_decline_clean:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[2].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)
'''
for sink_id in IDs_steady_decline_messy:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[2].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=alpha)
'''
axs.flatten()[2].set_ylabel("Mass accrretion rate ($\mathrm{M}_\odot/yr$)")

for sink_id in IDs_late_accretion_clean:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[3].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=0.25)
    
'''
for sink_id in IDs_late_accretion_messy:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[3].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=0.25)
'''
    
for sink_id in IDs_Other:
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data['dm'].T[sink_id]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[sink_id] - particle_data['tflush'].T[sink_id])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = (particle_data['time'].T[sink_id]-particle_data['time'].T[sink_id][0]) * units['time_unit'].in_units('Myr')
    axs.flatten()[4].semilogy(time_arr, acc_val, linewidth=ln_wd, alpha=0.25)
    
axs.flatten()[4].set_xlabel("Time since star formation (Myr)")

plt.savefig("FF_accretion_histories.jpg", bbox_inches='tight', pad_inches=0.02, dpi=300)
print("remade future fellows plot")
