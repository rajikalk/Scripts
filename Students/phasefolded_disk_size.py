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

#Phasefold
end_ind = np.argmin(abs(metadata['91'].age-yt.YTQuantity(120, 'kyr')))
time = metadata['91'].age[:end_ind]
separation = metadata['91'].separation[:end_ind]
disk_secondary = metadata['91'].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds)):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
    
bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_sec_phasefolded_all.png')

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds[:10])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
                
#plt.plot(t_scaled, disk_orb)

bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_sec_phasefolded_all_10.png')

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds[:20])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
                
#plt.plot(t_scaled, disk_orb)

bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_sec_phasefolded_all_20.png')

end_ind = 1200
time = metadata['90'].age[:end_ind]
separation = metadata['90'].separation[:end_ind]
disk_secondary = metadata['90'].disk_size[:end_ind]
ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0]
apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0]

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds)):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
    
bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_primary_phasefolded_all.png')

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds[:10])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
    
bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_primary_phasefolded_all_10.png')

plt.clf()
t_bin = np.linspace(0, 1, 11)
bin_vals = [[], [], [], [], [], [], [], [], [], []]
for peri_ind in range(1, len(periastron_inds[:20])):
    t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
    disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
    for bin_it in range(1, len(t_bin)):
        for t_val_it in range(len(t_scaled)):
            if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                bin_vals[bin_it-1].append(disk_orb[t_val_it])
    
bin_medians = []
bin_centers = (t_bin[1:] + t_bin[:-1])/2
for bin_val in bin_vals:
    bin_medians.append(np.median(bin_val))

plt.plot(bin_centers, bin_medians)
plt.savefig('Disk_primary_phasefolded_all_20.png')
