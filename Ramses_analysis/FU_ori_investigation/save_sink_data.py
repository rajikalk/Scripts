import numpy as np
import pickle
import glob
from pyramses import rsink
import sys
import os
import yt
import yt.units
from yt.units import g, s, cm, Lsun
import gc

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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
gc.collect()

if args.update_pickle == 'True':
    updating = False
    if args.sink_number == None:
        last_n = int(sorted(glob.glob(path+"output*"))[-1].split("_")[-1])
        stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
        while os.path.exists(stars_output_file) == False:
            last_n = last_n - 1
            stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
        loaded_sink_data_last = rsink(last_n, datadir=path)
        del last_n
        sink_ind = np.argmin(loaded_sink_data_last['u'])
        del loaded_sink_data_last
    else:
        sink_ind = args.sink_number
    print("Sink_ind =", sink_ind)
    
    gc.collect()
    
    if os.path.isfile('particle_data_'+str(sink_ind)+'.pkl'):
        try:
            print("reading", 'particle_data_'+str(sink_ind)+'.pkl')
            file_open = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            print("Reading particle data")
            loaded_sink_data = rsink(datadir=path, all=True)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
            del loaded_sink_data
        except:
            os.system('cp '+save_dir+'particle_data_'+str(sink_ind)+'_tmp.pkl '+save_dir+'particle_data_'+str(sink_ind)+'.pkl ')
            print("reading", 'particle_data_'+str(sink_ind)+'.pkl')
            file_open = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            print("Reading particle data")
            loaded_sink_data = rsink(datadir=path, all=True)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
            del loaded_sink_data
    else:
        updating = True
        particle_data = {}
        particle_data.update({'particle_tag':[]})
        particle_data.update({'time':[]})
        particle_data.update({'mass':[]})
        particle_data.update({'mdot':[]})
        particle_data.update({'separation':[]})
        particle_data.update({'secondary_position':[]})
        particle_data.update({'secondary_velocity':[]})
        
        counter = 0
        sink_form_time = 0
        
    gc.collect()
        
    if updating == True:
        print("Reading particle data")
        loaded_sink_data = rsink(datadir=path, all=True)
        loaded_sink_data = loaded_sink_data[counter:]
        gc.collect()
        while counter < len(loaded_sink_data):
            sink_data = loaded_sink_data[0]
            loaded_sink_data = loaded_sink_data[1:]
            counter = counter + 1
            if np.remainder(counter, 1000) == 0:
                try:
                    os.remove(save_dir+'particle_data_'+str(sink_ind)+'.pkl')
                except:
                    print("pickle files doesn't exist yet")
                file = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'wb')
                pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
                file.close()
                os.system('cp '+save_dir+'particle_data_'+str(sink_ind)+'.pkl '+save_dir+'particle_data_'+str(sink_ind)+'_tmp.pkl ')
                print('read', counter, 'snapshots of sink particle data, and saved pickle')
            if len(sink_data['u']) > sink_ind:
                tags = [sink_ind]
                pos_second = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
                dx = sink_data['x']*units['length_unit'].in_units('au') - pos_second[0]
                dy = sink_data['y']*units['length_unit'].in_units('au') - pos_second[1]
                dz = sink_data['z']*units['length_unit'].in_units('au') - pos_second[2]
                sep = np.sqrt(dx**2 + dy**2 + dz**2)
                del dx, dy, dz
                gc.collect()
                nearest_sink = np.argsort(sep)[1]
                del sep
                gc.collect()
                tags.append(nearest_sink)
                
                for tag in tags:
                    if tag not in particle_data['particle_tag']:
                        particle_data['particle_tag'].append(tag)
                
                pos_prim = yt.YTArray(np.array([sink_data['x'][nearest_sink], sink_data['y'][nearest_sink], sink_data['z'][nearest_sink]])*units['length_unit'].in_units('au'), 'au')
                vel_second = yt.YTArray(np.array([sink_data['ux'][sink_ind], sink_data['uy'][sink_ind], sink_data['uz'][sink_ind]])*units['velocity_unit'].in_units('km/s'), 'km/s')
                particle_data['secondary_position'].append(pos_second)
                particle_data['secondary_velocity'].append(vel_second)
                del vel_second
                gc.collect()
                
                separation = np.sqrt(np.sum((pos_second - pos_prim)**2))
                particle_data['separation'].append(separation)
                del separation
                gc.collect()
                if sink_form_time == 0:
                    sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
                particle_data['time'].append(time_val)
                del time_val
                gc.collect()
                particle_data['mass'].append(yt.YTArray(sink_data['m'][np.array([sink_ind, nearest_sink])]*units['mass_unit'].in_units('msun'), 'msun'))
                
                d_mass = sink_data['dm'][np.array([sink_ind, nearest_sink])]*units['mass_unit'].in_units('msun')
                d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
                acc_val = d_mass/d_time
                del d_mass, d_time
                    gc.collect()
                acc_val[np.where(acc_val == 0)[0]]=1.e-12
                particle_data['mdot'].append(yt.YTArray(acc_val, 'msun/yr'))
        #write lastest pickle
        file = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'wb')
        pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
        file.close()
        os.system('cp '+save_dir+'particle_data_'+str(sink_ind)+'.pkl '+save_dir+'particle_data_'+str(sink_ind)+'_tmp.pkl ')
        print('read', counter, 'snapshots of sink particle data, and saved pickle')
             
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

if 'particle_data' not in locals():
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
    print("Sink_ind =", sink_ind)
    file_open = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()

f_acc = 0.5
radius = yt.YTQuantity(2.0, 'rsun')
#M_dot = accretion(sink_inds, global_ind)
#M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
m_dot = yt.YTArray(particle_data['mdot']).in_units('g/s')
mass = yt.YTArray(particle_data['mass']).in_units('g')
L_acc = f_acc * (yt.units.gravitational_constant_cgs * mass * m_dot)/radius.in_units('cm')
L_tot = L_acc.in_units('Lsun')

Mag = -2.5*np.log10(L_tot)

plt.clf()
for part_it in range(2):
    if part_it == 0:
        plt.plot(particle_data['time'], Mag.T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
    else:
        plt.plot(particle_data['time'], Mag.T[part_it], label="Nearest_sink")
plt.gca().invert_yaxis()
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Magnitude (M$_{bol}$)')
plt.legend()
plt.title('Sink no ' + str(sink_ind))
plt.savefig(str(sink_ind)+'_magnitude_vs_time_sink.png')

plt.clf()
for part_it in range(2):
    if part_it == 0:
        plt.plot(particle_data['time'], L_tot.T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
    else:
        plt.plot(particle_data['time'], L_tot.T[part_it], label="Nearest_sink")
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Luminosity (Lsun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig(str(sink_ind)+'_luminosity_vs_time_sink.png')


plt.clf()
for part_it in range(2):
    if part_it == 0:
        plt.plot(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
    else:
        plt.plot(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Nearest_sink")
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Accretion rate (Msun/yr)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig(str(sink_ind)+'_accretion_vs_time_sink.png')

plt.clf()
for part_it in range(2):
    if part_it == 0:
        plt.plot(particle_data['time'], np.array(particle_data['mass']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
    else:
        plt.plot(particle_data['time'], np.array(particle_data['mass']).T[part_it], label="Nearest_sink")
plt.xlabel('Time (yr)')
plt.xlim()
plt.legend()
plt.ylabel('Mass (Msun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig(str(sink_ind)+'_mass_vs_time_sink.png')

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
plt.savefig(str(sink_ind)+'_separation_vs_time_sink.png')

if sink_ind == 45:
    start_time = 3500
    start_ind = np.argmin(abs(np.array(particle_data['time']) - start_time))
    end_time = 8500
    end_ind = np.argmin(abs(np.array(particle_data['time']) - end_time))
    plt.clf()
    plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'][start_ind:end_ind])
    plt.xlabel('Time (yr)')
    plt.xlim()
    plt.ylim(bottom=1.e-7)
    plt.ylabel('Accretion rate (Msun/yr)')
    plt.title('Sink no ' + str(sink_ind))
    plt.tick_params(axis='both', which='major', right=True, direction='in')
    plt.tick_params(axis='both', which='minor', right=True, direction='in')
    plt.savefig(str(sink_ind)+'_accretion_vs_time_sink_trunc_3500.png')
