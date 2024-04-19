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

path_to_pickle='/lustre/astro/vitot/ds_results/parameter_ds0.8_0.17.pkl'

information=open(path_to_pickle,"rb")
metadata=pickle.load(information,encoding='bytes')
information.close()
print(metadata.keys())
print(vars(metadata['90']))
print(len(metadata['90'].separation ))
print(len(metadata['90'].disk_size ))
print(len(metadata['91'].separation ))
print(len(metadata['91'].disk_size ))


orb1 = 0
orb2 = 1
for key in metadata.keys():
    if key == '13':
        pass
    else:
        ds_left = (metadata[key].separation[1:-1].value - metadata[key].separation[:-2].value)/(metadata[key].age[1:-1].value - metadata[key].age[:-2].value)
        ds_right = (metadata[key].separation[2:].value - metadata[key].separation[1:-1].value)/(metadata[key].age[2:].value - metadata[key].age[1:-1].value)
        metadata[key].periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0] + 1 
        metadata[key].apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0] + 1
        ds = metadata[key].disk_size[metadata[key].periastron_inds[orb1]:metadata[key].periastron_inds[orb2]]
        print(' SINK ', key)
        #print(max(ds),min(ds),np.mean(ds),  np.std(ds), np.median(ds))
        print('------')
        
for sec, prim in [('91','90'),('48','49'),('165','164')]:
    
    l_t1 = metadata[prim].periastron_inds[orb1]
    h_t1 = metadata[prim].periastron_inds[orb2]


    l_t2 = metadata[sec].periastron_inds[orb1]
    h_t2 = metadata[sec].periastron_inds[orb2]
    offset = np.abs(metadata[prim].age[l_t1]-metadata[sec].age[l_t2])
    if sec == '48':
        metadata[sec].age=metadata[sec].age-offset
    else:
        metadata[sec].age=metadata[sec].age+offset


#==================================================
two_col_width = 7.20472 #inches 
page_height = 10.62472 #inches
font_size = 8
c=1
plt.rcParams.update({'font.size': font_size})

# plt.rcParams['axes.labelsize'] = font_size
# plt.rcParams['mathtext.fontset'] = 'stixsans'
# plt.rcParams['mathtext.it'] = 'Arial:italic'
# plt.rcParams['mathtext.rm'] = 'Arial'
# plt.rcParams['mathtext.bf'] = 'Arial:bold'
# plt.rcParams['mathtext.it'] = 'Arial:italic'
# plt.rcParams['mathtext.rm'] = 'Arial'
# plt.rcParams['mathtext.sf'] = 'Arial'
# plt.rcParams['mathtext.default'] = 'regular'
# plt.rcParams['font.sans-serif'] = 'Arial'
# #plt.rcParams['font.family'] = 'sans-serif'

two_col_width = 7.20472 #inches 
single_col_width = 3.50394 #inches 
page_height = 10.62472 #inches


#Phasefold
min_orb = 1
max_orb = 11
n_bins = 20
plt.clf()
fig,axs = plt.subplots(2,3,tight_layout=True,figsize=(two_col_width,page_height/3))

ax=axs.flatten()
for i,sink_inds in enumerate([('91','90'),('48','49'),('165','164')]):
    end_ind = np.argmin(abs(metadata[sink_inds[0]].age.value-yt.YTQuantity(120, 'kyr').value))
    if i + 3 == 4:
        end_ind=-1
    end_ind=-1

    time = metadata[sink_inds[0]].age[:end_ind].value
    separation = metadata[sink_inds[0]].separation[:end_ind].value
    disk_secondary = metadata[sink_inds[0]].disk_size[:end_ind].value
    ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
    ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
    periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0] + 1 
    apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0] + 1

    
    t_bin = np.linspace(0, 1, n_bins +1)
    bin_mean_vals = []
    bin_median_vals = []
    bins_all = []
    for bin_it in range(1, len(t_bin)):
        bin_mean_vals.append([])
        bin_median_vals.append([])
        bins_all.append([])
    for peri_ind in range(min_orb, len(periastron_inds[min_orb:max_orb])):
        
        t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
        t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
        disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
        disk_orb[disk_orb==np.nan]=0
        for bin_it in range(1, len(t_bin)):
            bin_sub_set = []
            for t_val_it in range(len(t_scaled)):
                if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                    bin_sub_set.append(disk_orb[t_val_it])
            if np.isnan(np.median(bin_sub_set)) == False:
                bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
                bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
                bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set

    bin_medians = []
    bin_means = []
    bin_errs = []
    
    bin_centers = (t_bin[1:] + t_bin[:-1])/2
    for bin_val in bin_median_vals:
        median = np.median(bin_val)
        mean = np.mean(bin_val)
        std = np.std(bin_val)
        err = [median-(mean-std), (mean+std)-median]
        bin_means.append(mean)
        bin_medians.append(median)
        bin_errs.append(err)
            
    if sink_inds ==('91','90'):
        print('bin_centers',bin_centers) 
        print('bin_means',bin_means) 
        print('error',np.array(bin_errs).T)
    ax[i+3].errorbar(bin_centers, bin_means, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='secondary',color='blue')
    ax[i+3].set_ylim(bottom=0)
    #plt.plot(bin_centers, bin_medians)

    end_ind = 1200
    if i == 1:
        end_ind = -1
    end_ind = -1
    time = metadata[sink_inds[1]].age[:end_ind].value
    separation = metadata[sink_inds[1]].separation[:end_ind].value
    disk_secondary = metadata[sink_inds[1]].disk_size[:end_ind].value
    ds_left = (separation[1:-1] - separation[:-2])/(time[1:-1] - time[:-2])
    ds_right = (separation[2:] - separation[1:-1])/(time[2:] - time[1:-1])
    periastron_inds = np.argwhere((ds_left<0)&(ds_right>0)).T[0] + 1 
    apastron_inds = np.argwhere((ds_left>0)&(ds_right<0)).T[0] + 1

    bin_mean_vals = []
    bin_median_vals = []
    bins_all = []
    all_bin_ds=[]

    for bin_it in range(1, len(t_bin)):
        bin_mean_vals.append([])
        bin_median_vals.append([])
        bins_all.append([])
        all_bin_ds.append([])
    for peri_ind in range(min_orb, len(periastron_inds[min_orb:max_orb])):
        t_orb = time[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
        t_scaled = (t_orb - t_orb[0])/((t_orb - t_orb[0])[-1])
        disk_orb = disk_secondary[periastron_inds[peri_ind-1]:periastron_inds[peri_ind]]
        disk_orb[disk_orb==np.nan]=0

        for bin_it in range(1, len(t_bin)):
            bin_sub_set = []
            for t_val_it in range(len(t_scaled)):
                if t_scaled[t_val_it] >= t_bin[bin_it-1] and t_scaled[t_val_it] < t_bin[bin_it]:
                    bin_sub_set.append(disk_orb[t_val_it])
            all_bin_ds.append(bin_sub_set)
            if np.isnan(np.median(bin_sub_set)) == False:
                #excluding nans

                bin_mean_vals[bin_it-1].append(np.mean(bin_sub_set))
                bin_median_vals[bin_it-1].append(np.median(bin_sub_set))
                bins_all[bin_it-1] = bins_all[bin_it-1] + bin_sub_set
    #print(bin_sub_set)
    #print(np.shape(bin_sub_set))
    if sink_inds ==('91','90'):
        print('bin_median_vals',bin_median_vals)

    bin_medians = []
    bin_means = []

    bin_errs = []
    bin_centers = (t_bin[1:] + t_bin[:-1])/2
    for bin_val in bin_median_vals:
        median = np.median(bin_val)
        mean = np.mean(bin_val)
        std = np.std(bin_val)
        err = [median-(mean-std), (mean+std)-median]
        bin_means.append(mean)
        bin_medians.append(median)
        bin_errs.append(err)
        #print('Sink IDSS',sink_inds)

    ax[i].errorbar(bin_centers, bin_means, yerr=np.array(bin_errs).T, drawstyle='steps-mid', alpha=0.5, label='primary',color='orangered')
    if i == 0:
        anchored_text = AnchoredText('B1', 
                        frameon=False, borderpad=0, pad=0.1, 
                        loc=2, bbox_to_anchor=[0.9,0.96], 
                        bbox_transform=ax[0].transAxes,
                        prop={'fontsize':font_size+c,'fontfamily':'sans-serif','fontweight':'bold'})
        ax[0].add_artist(anchored_text)

        anchored_text = AnchoredText('B2', 
                        frameon=False, borderpad=0, pad=0.1, 
                        loc=2, bbox_to_anchor=[0.9,0.96], 
                        bbox_transform=ax[1].transAxes,
                        prop={'fontsize':font_size+c,'fontfamily':'sans-serif','fontweight':'bold'})
        ax[1].add_artist(anchored_text)

        anchored_text = AnchoredText('B3', 
                        frameon=False, borderpad=0, pad=0.1, 
                        loc=2, bbox_to_anchor=[0.9,0.96], 
                        bbox_transform=ax[2].transAxes,
                        prop={'fontsize':font_size+c,'fontfamily':'sans-serif','fontweight':'bold'})
        ax[2].add_artist(anchored_text)
    ax[0].legend(loc='upper left')
    ax[3].legend(loc='upper left')
    ax[i].set_ylim(bottom=0)
    ax[i+3].set_xlabel(r'phase $(\phi) $')
    ax[0].set_ylabel(r'd [AU]')
    ax[3].set_ylabel(r'd [AU]')
    #plt.savefig('/lustre/astro/vitot/20_MEANS_figure7_1_11.pdf',dpi=350)