import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pyramses import rsink
import sys
import os
import yt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-pickle", "--sink_pickle", type=str, default='particle_data.pkl')
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)
    
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})


file_open = open(save_dir+args.sink_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

f_acc = 0.5
radius = yt.YTQuantity(2.0, 'rsun')
m_dot = yt.YTArray(particle_data['mdot']).in_units('g/s')
mass = yt.YTArray(particle_data['mass']).in_units('g')
L_acc = f_acc * (yt.units.gravitational_constant * mass * m_dot)/radius.in_units('cm')
L_tot = L_acc.in_units('Lsun')

Mag = -2.5*np.log10(L_tot)
#Mag = np.nan_to_num(Mag)

separation = np.array(particle_data['separation'])
time = np.array(particle_data['time'])

ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
periastron_inds = periastron_inds + 1
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]
apastron_inds = apastron_inds + 1

#Make diagnostic plot
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width))
plt.semilogy(particle_data['time'], particle_data['separation'])
for peri_ind in periastron_inds:
    plt.axvline(x=np.array(particle_data['time'])[peri_ind])
plt.xlim([3500, 8700])
plt.savefig('periastrons.pdf')

#Find pre-peri time
pre_time = 10
pre_inds = []
end_inds = []
for peri_ind in range(len(periastron_inds)):
    target_time_start = time[periastron_inds[peri_ind]] - pre_time
    target_time_end = time[periastron_inds[peri_ind]] + 300
    try:
        if target_time_end > time[periastron_inds[peri_ind+1]]:
            target_time_end = time[periastron_inds[peri_ind+1]]
    except:
        pass
    pre_ind = np.argmin(abs(time - target_time_start))
    end_ind = np.argmin(abs(time - target_time_end))
    pre_inds.append(pre_ind)
    end_inds.append(end_ind)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, single_col_width))
ylim = [0, 7]
sig_thres = 3
plot_orbits = [1, 3, 4, 5]
for orb_it in plot_orbits: #range(5):
    time_orb = time[periastron_inds[orb_it]:periastron_inds[orb_it+1]]
    phase_orb = (time_orb - time_orb[0])/(time_orb - time_orb[0])[-1]
    acc_orb = m_dot[periastron_inds[orb_it]:periastron_inds[orb_it+1]].T[1].in_units('msun/yr')
    scaled_acc = acc_orb/np.max(acc_orb[:int(len(acc_orb)/2)])

    '''
    time_orb = time[pre_inds[orb_it-1]: end_inds[orb_it]] - time[periastron_inds[orb_it-1]]
    Mag_orb = m_dot[pre_inds[orb_it-1]: end_inds[orb_it]].in_units('msun/yr')
    Mag_orb_bounds = Mag_orb
    Mag_orb_bounds[np.where(Mag_orb_bounds==np.inf)] = np.nan
    mag_low = np.nanmax(Mag_orb_bounds)
    mag_std = np.nanstd(Mag_orb_bounds)
    mag_median = np.nanmedian(Mag_orb_bounds)
    mag_mean = np.nanmean(Mag_orb_bounds)
    #mag_sig = (mag_low - mag_median)/mag_std
    #mag_sig = (mag_low - mag_mean)/mag_std
    #if mag_sig > sig_thres:
    #    if np.nanmin(Mag_orb_bounds) < np.min(ylim):
    #        ylim = [np.nanmin(Mag_orb_bounds), ylim[1]]
    #    if np.nanmax(Mag_orb_bounds) > np.max(ylim):
    #        ylim = [ylim[0], np.nanmax(Mag_orb_bounds)]
        
    Mag_orb[np.where(np.isnan(Mag_orb) == True)] = np.inf
    '''
    plt.semilogy(phase_orb, scaled_acc, label="Orbit "+str(orb_it+1), color=colors[orb_it])
    #plt.semilogy(phase_orb, acc_orb.T[1], color=colors[orb_it])
    
plt.xlabel("Phase")
plt.ylabel("Mdot/Mdot_max")
plt.ylim([1.e-2, 1])
#plt.gca().invert_yaxis()
plt.xlim([0, 0.5])
plt.legend(loc='best')
plt.savefig('accretion_over_phase.png', bbox_inches='tight', dpi=300, pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, single_col_width))
ylim = [0, 7]
sig_thres = 3
#plot_orbits = [1, 3, 4, 5]
for orb_it in range(5): #plot_orbits: #range(5):
    time_orb = time[pre_inds[orb_it]:end_inds[orb_it]]
    time_orb = time_orb - time_orb[0]
    acc_orb = m_dot[pre_inds[orb_it]:end_inds[orb_it]].T[1].in_units('msun/yr')
    scaled_acc = acc_orb/np.max(acc_orb)

    '''
    time_orb = time[pre_inds[orb_it-1]: end_inds[orb_it]] - time[periastron_inds[orb_it-1]]
    Mag_orb = m_dot[pre_inds[orb_it-1]: end_inds[orb_it]].in_units('msun/yr')
    Mag_orb_bounds = Mag_orb
    Mag_orb_bounds[np.where(Mag_orb_bounds==np.inf)] = np.nan
    mag_low = np.nanmax(Mag_orb_bounds)
    mag_std = np.nanstd(Mag_orb_bounds)
    mag_median = np.nanmedian(Mag_orb_bounds)
    mag_mean = np.nanmean(Mag_orb_bounds)
    #mag_sig = (mag_low - mag_median)/mag_std
    #mag_sig = (mag_low - mag_mean)/mag_std
    #if mag_sig > sig_thres:
    #    if np.nanmin(Mag_orb_bounds) < np.min(ylim):
    #        ylim = [np.nanmin(Mag_orb_bounds), ylim[1]]
    #    if np.nanmax(Mag_orb_bounds) > np.max(ylim):
    #        ylim = [ylim[0], np.nanmax(Mag_orb_bounds)]
        
    Mag_orb[np.where(np.isnan(Mag_orb) == True)] = np.inf
    '''
    plt.semilogy(time_orb, scaled_acc, label="Orbit "+str(orb_it+1), color=colors[orb_it])
    #plt.semilogy(phase_orb, acc_orb.T[1], color=colors[orb_it])
    
plt.xlabel("Phase")
plt.ylabel("Mdot/Mdot_max")
plt.ylim([1.e-3, 1])
#plt.gca().invert_yaxis()
plt.xlim([-10, 300])
plt.legend(loc='best')
plt.savefig('accretion_over_time.png', bbox_inches='tight', dpi=300, pad_inches=0.02)
    
plot_orbits = [1, 3, 4, 5]
burst_time = [[5670, 5700], [7320, 7350], [7855, 7900], [8275, 8295]]
suppression_event = [[5575, 5700], [7290, 7350], [7850, 7900], [8265, 8295]]
for orb_it in range(len(plot_orbits)):
    plt.clf()
    time_orb = time[pre_inds[plot_orbits[orb_it]]:end_inds[plot_orbits[orb_it]]]
    #time_orb = time_orb - time_orb[0]
    acc_orb = m_dot[pre_inds[plot_orbits[orb_it]]:end_inds[plot_orbits[orb_it]]].T[1].in_units('msun/yr')
    scaled_acc = acc_orb/np.max(acc_orb)
    plt.semilogy(time_orb, scaled_acc)
    plt.axvspan(burst_time[orb_it][0], burst_time[orb_it][1], alpha=0.5)
    plt.axvspan(suppression_event[orb_it][0], suppression_event[orb_it][1], alpha=0.5)
    plt.xlabel("time since primary formation (yr)")
    plt.ylabel("scaled accretion")
    plt.xlim([time_orb[0], time_orb[-1]])
    plt.savefig('orbit_'+str(plot_orbits[orb_it]+1)+'.png', bbox_inches='tight', dpi=300, pad_inches=0.02)
