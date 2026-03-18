import yt
#yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
from mpi4py.MPI import COMM_WORLD as CW
import pickle
from pyramses import rsink
import yt.units
from yt.units import g, s, cm, Lsun
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

zoom_directory = "/home/100/rlk100/rlk/RAMSES/Zoom-in/"
Sim_dirs = [x[0] for x in os.walk(zoom_directory)]
Cleaned_dirs = []
for sim_dir in Sim_dirs:
    if sim_dir[-4:] == 'data':
        Cleaned_dirs.append(sim_dir)
del Sim_dirs

save_dir_top = "./"

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override

rit = -1
for path in Cleaned_dirs:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        sink_ind = int(path.split('Sink_')[-1].split('/')[0])
        pickle_name = 'particle_data_'+str(sink_ind)
        if "Event_restart" in path:
            save_dir = save_dir_top + 'Starting_from_event/'
        else:
            save_dir = save_dir_top
        if "Level" in path:
            lvl_int = path.split('Level_')[-1].split('/')[0]
            pickle_name = pickle_name + "_" + lvl_int
        pickle_name = pickle_name + ".pkl"
        print("Updating pickle for sim", path, "on rank", rank)

        updating = False
        
        if os.path.isfile(save_dir + pickle_name):
            try:
                print("reading", save_dir + pickle_name)
                file_open = open(save_dir + pickle_name, 'rb')
                particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
                file_open.close()
                counter = int(counter)
                print("Reading particle data")
                loaded_sink_data = rsink(datadir=path, all=True)
                if counter < len(loaded_sink_data):
                    updating = True
                    print('pickle data is not up to date! Updating...')
                #del loaded_sink_data
            except:
                os.system('cp '+save_dir + pickle_name.split('.pkl')[0]+'_tmp.pkl '+save_dir + pickle_name)
                print("reading", save_dir + pickle_name)
                file_open = open(save_dir + pickle_name, 'rb')
                particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
                file_open.close()
                counter = int(counter)
                print("Reading particle data")
                loaded_sink_data = rsink(datadir=path, all=True)
                if counter < len(loaded_sink_data):
                    updating = True
                    print('pickle data is not up to date! Updating...')
                #del loaded_sink_data
        else:
            updating = True
            particle_data = {}
            particle_data.update({'particle_tag':[]})
            particle_data.update({'time':[]})
            particle_data.update({'mass':[]})
            particle_data.update({'mdot':[]})
            particle_data.update({'separation':[]})
            loaded_sink_data = rsink(datadir=path, all=True)
            #particle_data.update({'secondary_position':[]})
            #particle_data.update({'secondary_velocity':[]})
            
            counter = 0
            sink_form_time = 0
            
        if updating == True:
            loaded_sink_data = loaded_sink_data[counter:]

            print("starting to update pickles, current counter=", counter)
            while len(loaded_sink_data)>0:
                sink_data = loaded_sink_data[0]
                loaded_sink_data = loaded_sink_data[1:]
                counter = counter + 1
                if np.remainder(counter, 100) == 0:
                    try:
                        os.remove(save_dir + pickle_name)
                    except:
                        print("pickle files doesn't exist yet")
                    file = open(save_dir + pickle_name, 'wb')
                    pickle.dump((particle_data, counter, sink_ind, sink_form_time), file, protocol=pickle.HIGHEST_PROTOCOL)
                    file.close()
                    os.system('cp '+save_dir + pickle_name + ' ' + save_dir+pickle_name.split('.pkl')[0]+'_tmp.pkl')
                    print('read', counter, 'snapshots of sink particle data, and saved pickle')
                    
                    #make progress plots
                    if len(particle_data['mass'])>0:
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
                        plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_magnitude_vs_time_sink.png')

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
                        plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_luminosity_vs_time_sink.png')


                        plt.clf()
                        for part_it in range(2):
                            if part_it == 0:
                                plt.semilogy(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
                            else:
                                plt.semilogy(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Nearest_sink")
                        plt.xlabel('Time (yr)')
                        plt.xlim()
                        plt.legend()
                        plt.ylabel('Accretion rate (Msun/yr)')
                        plt.title('Sink no ' + str(sink_ind))
                        plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_accretion_vs_time_sink.png')

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
                        plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_mass_vs_time_sink.png')

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
                        plt.title('Sink no ' + str(sink_ind) + " with companion tags " + str(particle_data['particle_tag'][1:]))
                        plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_separation_vs_time_sink.png')

                if len(sink_data['u']) > sink_ind:
                    tags = [sink_ind]
                    pos_second = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
                    dx = sink_data['x']*units['length_unit'].in_units('au') - pos_second[0]
                    dy = sink_data['y']*units['length_unit'].in_units('au') - pos_second[1]
                    dz = sink_data['z']*units['length_unit'].in_units('au') - pos_second[2]
                    sep = np.sqrt(dx**2 + dy**2 + dz**2)
                    #del dx, dy, dz
                    nearest_sink = np.argsort(sep)[1]
                    #del sep
                    tags.append(nearest_sink)
                    
                    for tag in tags:
                        if tag not in particle_data['particle_tag']:
                            particle_data['particle_tag'].append(tag)
                    
                    pos_prim = yt.YTArray(np.array([sink_data['x'][nearest_sink], sink_data['y'][nearest_sink], sink_data['z'][nearest_sink]])*units['length_unit'].in_units('au'), 'au')
                    #vel_second = yt.YTArray(np.array([sink_data['ux'][sink_ind], sink_data['uy'][sink_ind], sink_data['uz'][sink_ind]])*units['velocity_unit'].in_units('km/s'), 'km/s')
                    #particle_data['secondary_position'].append(pos_second)
                    #particle_data['secondary_velocity'].append(vel_second)
                    #del vel_second
                    
                    separation = np.sqrt(np.sum((pos_second - pos_prim)**2))
                    particle_data['separation'].append(separation)
                    #del separation

                    if sink_form_time == 0:
                        sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                    time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
                    particle_data['time'].append(time_val)
                    #del time_val

                    particle_data['mass'].append(yt.YTArray(sink_data['m'][np.array([sink_ind, nearest_sink])]*units['mass_unit'].in_units('msun'), 'msun'))
                    
                    d_mass = sink_data['dm'][np.array([sink_ind, nearest_sink])]*units['mass_unit'].in_units('msun')
                    d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
                    acc_val = d_mass/d_time
                    #del d_mass, d_time
                    
                    acc_val[np.where(acc_val == 0)[0]]=1.e-12
                    particle_data['mdot'].append(yt.YTArray(acc_val, 'msun/yr'))
                    #print('read', counter)
            #write lastest pickle
            file = open(save_dir + pickle_name, 'wb')
            pickle.dump((particle_data, counter, sink_ind, sink_form_time), file, protocol=pickle.HIGHEST_PROTOCOL)
            file.close()
            os.system('cp '+save_dir + pickle_name + ' '+save_dir+pickle_name.split('.pkl')[0]+'_tmp.pkl')
            print('read', counter, 'snapshots of sink particle data, and saved pickle')

        if 'particle_data' not in locals():
            file_open = open(save_dir + pickle_name, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()

        if len(particle_data['mass'])>0:
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
            plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_magnitude_vs_time_sink.png')

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
            plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_luminosity_vs_time_sink.png')


            plt.clf()
            for part_it in range(2):
                if part_it == 0:
                    plt.semilogy(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Sink "+str(particle_data['particle_tag'][part_it]))
                else:
                    plt.semilogy(particle_data['time'], np.array(particle_data['mdot']).T[part_it], label="Nearest_sink")
            plt.xlabel('Time (yr)')
            plt.xlim()
            plt.legend()
            plt.ylabel('Accretion rate (Msun/yr)')
            plt.title('Sink no ' + str(sink_ind))
            plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_accretion_vs_time_sink.png')

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
            plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_mass_vs_time_sink.png')

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
            plt.title('Sink no ' + str(sink_ind) + " with companion tags " + str(particle_data['particle_tag'][1:]))
            plt.savefig(save_dir +pickle_name.split('.pkl')[0].split('data_')[-1]+'_separation_vs_time_sink.png')

