import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib.pyplot as plt
import pickle

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

#Plot sink histories for classification
for sink_id in range(len(particle_data['dm'][0])):
    zero_inds = np.where(particle_data['time'].T[sink_id] == 0)[0]
    particle_data['time'].T[sink_id][zero_inds] = np.nan
    particle_data['dm'].T[sink_id][zero_inds] = np.nan
    d_mass = particle_data[sink_ind][10]*units['mass_unit'].in_units('msun')
    d_time = (particle_data['time'].T[0] - particle_data['tflush'].T[0])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    time_arr = particle_data['time'].T[sink_id] * scale_t.in_units('yr')
    plt.clf()
    plt.semilogy(time_arr, acc_val)
    plt.xlabel("time")
    plt.ylabel("Accretion rate")
    plt.ylim([1.e-9, 5.e-5])
    plt.title("Sink ID = "+str(sink_id))
    plt.savefig("Accretion_history_"+str(sink_id)+".png")
    print("saved accretion history for sink", sink_id)
    
    
IDs_continuous = [2]
IDs_suddenly_shut_off = [0, 1,]
IDs_steady_decline = [3, 4, 5, 6]
IDs_late_accretion = []
IDs_Other = []
