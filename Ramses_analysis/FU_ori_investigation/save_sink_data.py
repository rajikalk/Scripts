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
            particle_data, counter, sink_ind = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            #loaded_sink_data = loaded_sink_data[counter:]
            import pdb
            pdb.set_trace()
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating')
        except:
            os.system('cp '+save_dir+'particle_data_tmp.pkl '+save_dir+'particle_data.pkl ')
            file_open = open(save_dir+'particle_data.pkl', 'rb')
            particle_data, counter, sink_ind = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            import pdb
            pdb.set_trace()
            #loaded_sink_data = loaded_sink_data[counter:]
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating')
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
        
        counter = 0
        sink_form_time = 0
        
    if updating == True:
        for sink_data in loaded_sink_data:
            counter = counter + 1
            if np.remainder(counter, 1000) == 0:
                try:
                    os.remove(save_dir+'particle_data.pkl')
                except:
                    print("pickle files doesn't exist yet")
                file = open(save_dir+'particle_data.pkl', 'wb')
                pickle.dump((particle_data, counter, sink_ind), file)
                file.close()
                os.system('cp '+save_dir+'particle_data.pkl '+save_dir+'particle_data_tmp.pkl ')
                print('read', counter, 'snapshots of sink particle data, and saved pickle')
            if len(sink_data['u']) > sink_ind:
                particle_data['particle_tag'].append(sink_ind)
                if sink_form_time == 0:
                    sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
                #if len(particle_tags) == 1:
                particle_data['time'].append(time_val)
                particle_data['mass'].append(yt.YTArray(sink_data['m'][sink_ind]*units['mass_unit'].in_units('msun'), 'msun'))
                
                d_mass = sink_data['dm'][sink_ind]*units['mass_unit'].in_units('msun')
                d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
                particle_data['mdot'].append(yt.YTArray(d_mass/d_time, 'msun/yr'))
                

f_acc = 0.5
radius = yt.YTQuantity(2.0, 'rsun')
#M_dot = accretion(sink_inds, global_ind)
#M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
m_dot = yt.YTArray(particle_data['mdot']).in_units('g/s')
mass = yt.YTArray(particle_data['mass']).in_units('g')
L_acc = f_acc * (yt.units.G * mass * m_dot)/radius.in_units('cm')
L_tot = L_acc.in_units('Lsun')

plt.clf()
plt.semilogy(particle_data['time'], L_tot)
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Luminosity (Lsun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('luminosity_vs_time_sink_'+str(sink_ind)+'.png')


plt.clf()
plt.semilogy(particle_data['time'], particle_data['mdot'])
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Accretion rate (Msun/yr)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('accretion_vs_time_sink_'+str(sink_ind)+'.png')

plt.clf()
plt.plot(particle_data['time'], particle_data['mass'])
plt.xlabel('Time (yr)')
plt.xlim()
plt.ylabel('Mass (Msun)')
plt.title('Sink no ' + str(sink_ind))
plt.savefig('mass_vs_time_sink_'+str(sink_ind)+'.png')
