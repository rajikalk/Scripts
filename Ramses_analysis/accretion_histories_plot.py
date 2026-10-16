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
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width*1))#, sharey=True)

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
        axs.semilogy(time_arr, acc_val, linewidth=0.25, alpha=0.2)
if plot_individual == False:
    axs.set_xlabel("Simulation time (Myr)")
    axs.set_ylabel("Accretion rate ($\mathrm{M}_\odot/\mathrm{yr}$)")
    if sink_id == 0:
        axs.set_xlim([time_arr[0], time_arr[-1]])
    axs.set_ylim([1.e-9, 5.e-5])
    plt.savefig("All_accretion_histories.pdf", bbox_inches='tight', pad_inches=0.02)
    print("saved all accretion histories")

        
    
    
IDs_continuous = [2]
IDs_suddenly_shut_off = [0, 1,]
IDs_steady_decline = [3, 4, 5, 6]
IDs_late_accretion = []
IDs_Other = []
