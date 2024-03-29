import yt
import matplotlib.pyplot as plt
import pickle
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as patches
class binary_timeseries:
  def __init__(self,roche,separation,age,disk_size,disk_mass,mass_in_100,mass_in_5000,output,accretion,sink_mass):
    self.age = age
    self.disk_size = disk_size
    self.roche = roche
    self.separation = separation
    self.disk_mass = disk_mass
    self.mass_in_100 = mass_in_100
    self.mass_in_5000 = mass_in_5000
    self.output = output
    self.accretion = accretion
    self.sink_mass = sink_mass
    self.units =  {'distances': 'au' , 'masses':'Msun','age':'kyr','accretion':'Lsun'}

class single_timeseries:
  def __init__(self, age,disk_size,disk_mass,mass_in_100,mass_in_5000,output,accretion,sink_mass):
    self.age = age
    self.disk_size = disk_size
    self.disk_mass = disk_mass
    self.mass_in_100 = mass_in_100
    self.mass_in_5000 = mass_in_5000
    self.outout = output
    self.accretion = accretion
    self.sink_mass = sink_mass
    self.units =  {'distances': 'au' , 'masses':'Msun','age':'kyr','accretion':'Lsun'}

path_to_pickle='/lustre/astro/vitot/viscosity_analysis/collected_metadata/HR_collected_metadata.pkl'

information=open(path_to_pickle,"rb")
metadata=pickle.load(information,encoding='bytes')
information.close()
print(metadata.keys())
print(vars(metadata['90']))
print(len(metadata['90'].separation ))
print(len(metadata['90'].disk_size ))
print(len(metadata['91'].separation ))
print(len(metadata['91'].disk_size ))

N_orb = 10
sink_inds = ['91', '90']

#Phasefold
end_ind = np.argmin(abs(metadata[sink_inds[0]].age-yt.YTQuantity(120, 'kyr')))
time = metadata[sink_inds[0]].age[:end_ind]
separation = metadata[sink_inds[0]].separation[:end_ind]
disk_secondary = metadata[sink_inds[0]].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

plt.clf()
plt.plot(time, separation)
for peri_ind in periastron_inds:
    plt.axvline(x=time[peri_ind+1], color='b')
for ap_ind in apastron_inds:
    plt.axvline(x=time[ap_ind+1], color='r')
plt.xlim([60, 90])
plt.ylim([0, 200])
plt.savefig('separation.png')

plt.clf()
t_bin = np.linspace(0, 1, 21)
bin_mean_vals = []
bin_median_vals = []
bins_all = []
for bin_it in range(1, len(t_bin)):
    bin_mean_vals.append([])
    bin_median_vals.append([])
    bins_all.append([])
for peri_ind in range(1, len(periastron_inds[:N_orb])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        bin_sub_set = []
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] <= t_bin[bin_it]:
                bin_sub_set.append(disk_orb[t_val_it])
        if np.isnan(np.median(bin_sub_set)) == False:
            bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
            bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
            bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set
    
bin_medians = []
bin_errs = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_mean_vals:
    median = np.median(bin_val)
    mean = np.mean(bin_val)
    std = np.std(bin_val)
    err = [median-(mean-std), (mean+std)-median]
    bin_medians.append(median)
    bin_errs.append(err)

plt.errorbar(bin_centers, bin_medians, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='secondary')
#plt.plot(bin_centers, bin_medians)

end_ind = 1200
time = metadata[sink_inds[1]].age[:end_ind]
separation = metadata[sink_inds[1]].separation[:end_ind]
disk_secondary = metadata[sink_inds[1]].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

t_bin = np.linspace(0, 1, 21)
bin_mean_vals = []
bin_median_vals = []
bins_all = []
for bin_it in range(1, len(t_bin)):
    bin_mean_vals.append([])
    bin_median_vals.append([])
    bins_all.append([])
for peri_ind in range(1, len(periastron_inds[:N_orb])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        bin_sub_set = []
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] <= t_bin[bin_it]:
                bin_sub_set.append(disk_orb[t_val_it])
        if np.isnan(np.median(bin_sub_set)) == False:
            bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
            bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
            bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set
    
bin_medians = []
bin_errs = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_mean_vals:
    median = np.median(bin_val)
    mean = np.mean(bin_val)
    std = np.std(bin_val)
    err = [median-(mean-std), (mean+std)-median]
    bin_medians.append(median)
    bin_errs.append(err)


plt.errorbar(bin_centers, bin_medians, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='primary')
plt.ylim(bottom=0)
plt.legend()
plt.savefig('phasefolded_all_'+str(N_orb)+ '_' + sink_inds[0] +'.png')

#==================================================

#Phasefold
end_ind = np.argmin(abs(metadata[sink_inds[0]].age-yt.YTQuantity(120, 'kyr')))
time = metadata[sink_inds[0]].age[:end_ind]
separation = metadata[sink_inds[0]].separation[:end_ind]
disk_secondary = metadata[sink_inds[0]].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

plt.clf()
t_bin = np.linspace(0, 1, 21)
bin_mean_vals = []
bin_median_vals = []
bins_all = []
for bin_it in range(1, len(t_bin)):
    bin_mean_vals.append([])
    bin_median_vals.append([])
    bins_all.append([])
for peri_ind in range(1, len(periastron_inds[:N_orb])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        bin_sub_set = []
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] <= t_bin[bin_it]:
                bin_sub_set.append(disk_orb[t_val_it])
        if np.isnan(np.median(bin_sub_set)) == False:
            bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
            bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
            bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set
    
bin_medians = []
bin_errs = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_mean_vals:
    median = np.median(bin_val)
    mean = np.mean(bin_val)
    std = np.std(bin_val)
    err = [median-(mean-std), (mean+std)-median]
    bin_medians.append(median)
    bin_errs.append(err)

plt.errorbar(bin_centers, bin_medians, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='secondary')
#plt.plot(bin_centers, bin_medians)

end_ind = 1200
time = metadata[sink_inds[1]].age[:end_ind]
separation = metadata[sink_inds[1]].separation[:end_ind]
disk_secondary = metadata[sink_inds[1]].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

t_bin = np.linspace(0, 1, 21)
bin_mean_vals = []
bin_median_vals = []
bins_all = []
for bin_it in range(1, len(t_bin)):
    bin_mean_vals.append([])
    bin_median_vals.append([])
    bins_all.append([])
for peri_ind in range(1, len(periastron_inds[:N_orb])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        bin_sub_set = []
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] <= t_bin[bin_it]:
                bin_sub_set.append(disk_orb[t_val_it])
        if np.isnan(np.median(bin_sub_set)) == False:
            bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
            bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
            bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set
    
bin_medians = []
bin_errs = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_mean_vals:
    median = np.median(bin_val)
    mean = np.mean(bin_val)
    std = np.std(bin_val)
    err = [median-(mean-std), (mean+std)-median]
    bin_medians.append(median)
    bin_errs.append(err)


plt.errorbar(bin_centers, bin_medians, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='primary')
plt.ylim(bottom=0)
plt.legend()
plt.savefig('phasefolded_median_all_'+str(N_orb)+ '_' + sink_inds[0] + '.png')

