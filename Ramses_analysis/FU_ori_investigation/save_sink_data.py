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
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()

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

if args.update_pickle == 'True':
    print("Reading particle data")
    loaded_sink_data = rsink(datadir=path, all=True)
    updating = False
    
    if os.path.isfile('particle_data.pkl'):
        try:
            file_open = open(save_dir+'particle_data.pkl', 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
        except:
            os.system('cp '+save_dir+'particle_data_tmp.pkl '+save_dir+'particle_data.pkl ')
            file_open = open(save_dir+'particle_data.pkl', 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
    else:
        updating = True
        if args.sink_number == None:
            last_n = int(sorted(glob.glob(path+"output*"))[-1].split("_")[-1])
            stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
            while os.path.exists(stars_output_file) == False:
                last_n = last_n - 1
                stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
            loaded_sink_data_last = rsink(last_n, datadir=path)
            sink_ind = np.argmin(loaded_sink_data_last['u'])
        else:
            sink_ind = args.sink_number
            
        particle_data = {}
        particle_data.update({'particle_tag':[]})
        particle_data.update({'time':[]})
        particle_data.update({'mass':[]})
        particle_data.update({'mdot':[]})
        particle_data.update({'separation':[]})
        
        counter = 0
        sink_form_time = 0
        
    if updating == True:
        loaded_sink_data = loaded_sink_data[counter:]
        for sink_data in loaded_sink_data:
            counter = counter + 1
            if np.remainder(counter, 1000) == 0:
                try:
                    os.remove(save_dir+'particle_data.pkl')
                except:
                    print("pickle files doesn't exist yet")
                file = open(save_dir+'particle_data.pkl', 'wb')
                pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
                file.close()
                os.system('cp '+save_dir+'particle_data.pkl '+save_dir+'particle_data_tmp.pkl ')
                print('read', counter, 'snapshots of sink particle data, and saved pickle')
            if len(sink_data['u']) > sink_ind:
                tags = np.arange(len(sink_data['u']))[sink_ind-1:sink_ind+1]
                for tag in tags:
                    if tag not in particle_data['particle_tag']:
                        particle_data['particle_tag'].append(tag)
                pos_prim = yt.YTArray(np.array([sink_data['x'][sink_ind-1], sink_data['y'][sink_ind-1], sink_data['z'][sink_ind-1]])*units['length_unit'].in_units('au'), 'au')
                pos_second = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
                separation = np.sqrt(np.sum((pos_second - pos_prim)**2))
                particle_data['separation'].append(separation)
                if sink_form_time == 0:
                    sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
                #if len(particle_tags) == 1:
                particle_data['time'].append(time_val)
                particle_data['mass'].append(yt.YTArray(sink_data['m'][sink_ind-1:sink_ind+1]*units['mass_unit'].in_units('msun'), 'msun'))
                
                d_mass = sink_data['dm'][sink_ind-1:sink_ind+1]*units['mass_unit'].in_units('msun')
                d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
                acc_val = d_mass/d_time
                if len(acc_val) == 1:
                    import pdb
                    pdb.set_trace()
                acc_val[np.where(acc_val == 0)[0]]=1.e-12
                particle_data['mdot'].append(yt.YTArray(acc_val, 'msun/yr'))
        #write lastest pickle
        file = open(save_dir+'particle_data.pkl', 'wb')
        pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
        file.close()
        os.system('cp '+save_dir+'particle_data.pkl '+save_dir+'particle_data_tmp.pkl ')
        print('read', counter, 'snapshots of sink particle data, and saved pickle')
                

f_acc = 0.5
radius = yt.YTQuantity(2.0, 'rsun')
#M_dot = accretion(sink_inds, global_ind)
#M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
m_dot = yt.YTArray(particle_data['mdot']).in_units('g/s')
mass = yt.YTArray(particle_data['mass']).in_units('g')
L_acc = f_acc * (yt.units.G * mass * m_dot)/radius.in_units('cm')
L_tot = L_acc.in_units('Lsun')

Mag = -2.5*np.log10(L_tot)

plt.clf()
for part_it in range(len(particle_data['particle_tag'])):
    plt.plot(particle_data['time'], Mag.T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]+1))
plt.gca().invert_yaxis()
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Magnitude (M$_{bol}$)')
plt.legend()
plt.title('Sink no ' + str(sink_ind))
plt.savefig('magnitude_vs_time_sink_'+str(sink_ind)+'.png')

plt.clf()
for part_it in range(len(particle_data['particle_tag'])):
    plt.semilogy(particle_data['time'], L_tot.T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]+1))
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Luminosity (Lsun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('luminosity_vs_time_sink_'+str(sink_ind)+'.png')


plt.clf()
for part_it in range(len(particle_data['particle_tag'])):
    plt.semilogy(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]+1))
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Accretion rate (Msun/yr)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('accretion_vs_time_sink_'+str(sink_ind)+'.png')

plt.clf()
for part_it in range(len(particle_data['particle_tag'])):
    plt.plot(particle_data['time'], np.array(particle_data['mass']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]+1))
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Mass (Msun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('mass_vs_time_sink_'+str(sink_ind)+'.png')

L = yt.YTQuantity(4, 'pc')
curr_dir = os.getcwd()
if 'Level' not in curr_dir:
    refinement = 18
else:
    refinement = int(curr_dir.split('Level_')[-1][:2])
d_min = L.in_units('au')/(2**refinement)
plt.clf()
plt.semilogy(particle_data['time'], particle_data['separation'])
plt.axhline(y=8*d_min, color='r', linestyle=':')
plt.axhline(y=4*d_min, color='r', linestyle='--')
plt.axhline(y=d_min, color='r')
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Separation (AU)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('separation_vs_time_sink_'+str(sink_ind)+'.png')

if sink_ind == 45:
    start_time = 3500
    start_ind = np.argmin(abs(np.array(particle_data['time']) - start_time))
    plt.clf()
    plt.semilogy(particle_data['time'][start_ind:], particle_data['mdot'][start_ind:])
    plt.xlabel('Time (yr)')
    plt.xlim()
    plt.ylim(bottom=1.e-7)
    plt.ylabel('Accretion rate (Msun/yr)')
    plt.title('Sink no ' + str(sink_ind))
    plt.tick_params(axis='both', which='major', right=True, direction='in')
    plt.tick_params(axis='both', which='minor', right=True, direction='in')
    plt.savefig('accretion_vs_time_sink_'+str(sink_ind)+'_trunc_3500.png')
