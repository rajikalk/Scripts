#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import os
import yt
yt.enable_parallelism()
import gc
from mpi4py.MPI import COMM_WORLD as CW

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

parser = argparse.ArgumentParser()
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=10000, type=float)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------

#Start by loading pickel data and then deleting what we don't need
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/data/'
files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))

if os.path.exists('Kep_mass.pkl'):
    file_open = open('Kep_mass.pkl', 'rb')
    save_dict = pickle.load(file_open)
    file_open.close()
    if len(save_dict["Time"]) != len(files):
        files = files[len(save_dict["Time"]):]
elif os.path.exists('Kep_mass_0.pkl'):
    pickle_files = sorted(glob.glob("Kep_mass_*.pkl"))
    save_dict = {}
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        save_dict_r = pickle.load(file_open)
        file_open.close()
        for key in save_dict_r.keys():
            if key not in save_dict.keys():
                save_dict.update({key:save_dict_r[key]})
            else:
                save_dict[key] = np.append(save_dict[key], save_dict_r[key])
    del save_dict_r
    gc.collect()
    sorted_inds = np.argsort(save_dict["Time"])
    for key in save_dict.keys():
        save_dict[key] = save_dict[key][sorted_inds]
    if len(save_dict["Time"]) != len(files):
        files = files[len(save_dict["Time"]):]
else:
    save_dict = {}


sink_id = 45
sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()
gc.collect()

#Define radius bins
radius = yt.YTQuantity(args.measuring_sphere_radius, 'au')
radius_bins = np.logspace(0, np.log10(radius), 100)

if len(files)>0:
    for fn in yt.parallel_objects(files, njobs=1):
        proj_root_rank = int(rank/para_div)
        print('Reading file', fn, 'on rank', rank)
        sys.stdout.flush()
        ds = yt.load(fn, units_override=units_override)
        
        if len(ds.r["sink_particle_form_time"]) == 45:
            skip=True
        else:
            if np.isnan(sink_form_time):
                sink_form_time = ds.r["sink_particle_form_time"][sink_id]
                print("RANK", rank, "got sink formation time")
                sys.stdout.flush()
            skip = False
        if skip == False:
            time_val = ds.current_time.in_units('yr').value - sink_form_time.in_units('yr').value
            if "Time" not in save_dict.keys():
                save_dict.update({"Time":np.array([time_val])})
            else:
                save_dict["Time"] = np.append(save_dict["Time"],time_val)
            del time_val
            gc.collect()
            print("RANK", rank, "Got time stamp")
            sys.stdout.flush()
            
            #Get sink position
            sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
            sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
            sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
            sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
            
            dx_sinks = ds.r["gas", "sink_particle_posx"].in_units('au') - sink_pos[0].in_units('au')
            dy_sinks = ds.r["gas", "sink_particle_posy"].in_units('au') - sink_pos[1].in_units('au')
            dz_sinks = ds.r["gas", "sink_particle_posz"].in_units('au') - sink_pos[2].in_units('au')
            sink_separations = np.sqrt(dx_sinks**2 + dy_sinks**2 + dz_sinks**2)
            del sink_particle_posx, sink_particle_posy, sink_particle_posz, dx_sinks, dy_sinks, dz_sinks
            gc.collect()
            print("RANK", rank, "Got sink position")
            sys.stdout.flush()
            
            #Get inds in measuring sphere
            dx = ds.r['ramses', 'x'].in_units('au') - sink_pos[0].in_units('au')
            dy = ds.r['ramses', 'y'].in_units('au') - sink_pos[1].in_units('au')
            dz = ds.r['ramses', 'z'].in_units('au') - sink_pos[2].in_units('au')
            sep_vector_all = yt.YTArray([dx, dy, dz])
            del dx, dy, dz, sink_pos
            gc.collect()
            print("RANK", rank, "Got separation vectors")
            sys.stdout.flush()
            sep = np.sqrt(sep_vector_all[0]**2 + sep_vector_all[1]**2 + sep_vector_all[2]**2)
            
            sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
            sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
            sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
            sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
            del sink_particle_velx, sink_particle_vely, sink_particle_velz
            gc.collect()
            print("RANK", rank, "got particle velocity")
            
            #get separations to other sink particles
            
            #Get indices in measure sphere
            #Start iterating over Radial bins
            for radius_bit in range(1, len(radius_bins)):
                #Calculate enclosed mass:
                enclosed_inds = np.where(sep<=radius_bins[radius_bit])[0]
                enclosed_mass = np.sum(ds.r["gas", "mass"][enclosed_inds])
                enclosed_sinks = np.where(sink_separations<=radius_bins[radius_bit])[0]
                enclosed_sink_mass = np.sum(ds.r["gas", "sink_particle_mass"][enclosed_sinks])
                enclosed_mass = enclosed_mass + enclosed_sink_mass
                del enclosed_sinks, enclosed_sink_mass, enclosed_inds
                gc.collect()
                
                #Now get indices in sphere
                sphere_inds = np.where((sep>radius_bins[radius_bit-1])&(sep<=radius_bins[radius_bit]))[0]
                #get average radius in the bin
                shell_radius = np.mean(sep[sphere_inds])
                #calcualte gravitational potential energy
                E_grav = (yt.units.gravitational_constant_cgs*enclosed_mass*ds.r["gas", "mass"][sphere_inds])/shell_radius
                #calcualte kinetic energy
                sph_dvx = ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s') - sink_vel[0]
                sph_dvy = ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s') - sink_vel[1]
                sph_dvz = ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s') - sink_vel[2]
                sph_vel = yt.YTArray([sph_dvx, sph_dvy, sph_dvz])
                E_kin = 0.5 * ds.r["gas", "mass"][sphere_inds] * sph_vel**2
                
                import pdb
                pdb.set_trace()
                
                
            radii = sep[sphere_inds]
            del sep
            gc.collect()
            print("RANK", rank, "Got indices in measuring sphere")
            sys.stdout.flush()
            sep_vector = sep_vector_all.T[sphere_inds].T
            del sep_vector_all
            gc.collect()
            print("RANK", rank, "Got separation vectors")
            sys.stdout.flush()
            
            #Calcualte enclosed mass
            try:
                gas_mass = ds.r["gas", "mass"][sphere_inds]
            except:
                gas_mass = ds.r["gas", "Density"][sphere_inds].in_units('g/cm**3')*(ds.r["ramses", "dx"][sphere_inds].in_units('cm')**3)
            enclosed_mass = yt.YTArray(np.zeros(np.shape(radii)), "g")
            for radi_it in range(len(radii)):
                enc_inds = np.where(radii<=radii[radi_it])[0]
                enc_mass = np.sum(gas_mass[enc_inds])
                enclosed_mass[radi_it] = enc_mass
            gc.collect()
            print("RANK", rank, "calculated enclosed gas mass")
            sys.stdout.flush()
            enclosed_mass = enclosed_mass+sink_mass.in_units('g')
            del sink_mass
            gc.collect()
            print("RANK", rank, "got sink mass")
            sys.stdout.flush()
            keplerian_velocity = np.sqrt((yt.units.gravitational_constant_cgs*enclosed_mass)/radii).in_units('km/s')
            del enclosed_mass
            gc.collect()
            print("RANK", rank, "calculated keplerian mass")
            sys.stdout.flush()
            
            #Calculate bulk velocity of the sphere
            sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
            sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
            sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
            sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
            del sink_particle_velx, sink_particle_vely, sink_particle_velz
            gc.collect()
            print("RANK", rank, "got particle velocity")
            sys.stdout.flush()
            #print('Got particle velocity on rank', rank, ' for fn', ds)
            sys.stdout.flush()
            sph_dvx = ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s') - sink_vel[0]
            sph_dvy = ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s') - sink_vel[1]
            sph_dvz = ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s') - sink_vel[2]
            sph_vel = yt.YTArray([sph_dvx, sph_dvy, sph_dvz])
            del sph_dvx, sph_dvy, sph_dvz, sink_vel
            gc.collect()
            print("RANK", rank, "got gas mass in measuring sphere")
            sys.stdout.flush()
            
            #Calculate Tangential vel
            proj_factor = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()
            sph_speed = np.sqrt(sph_vel[0]**2 + sph_vel[1]**2 + sph_vel[2]**2)
            rel_kep_full = sph_speed/keplerian_velocity
            del sph_vel
            gc.collect()
            print("RANK", rank, "got gas speed and proj factor")
            sys.stdout.flush()

            proj_v_x = proj_factor*sep_vector.in_units('km')[0]
            proj_v_y = proj_factor*sep_vector.in_units('km')[1]
            proj_v_z = proj_factor*sep_vector.in_units('km')[2]
            del proj_factor, sep_vector
            gc.collect()
            print("RANK", rank, "calculated projected components")
            sys.stdout.flush()
            rad_speed = np.sqrt(proj_v_x**2 + proj_v_y**2 + proj_v_z**2)
            del proj_v_x, proj_v_y, proj_v_z
            gc.collect()
            print("RANK", rank, "calculated radial velocity")
            sys.stdout.flush()
            
            tang_vel = np.sqrt(sph_speed**2 - rad_speed**2)
            del sph_speed, rad_speed
            gc.collect()
            print("RANK", rank, "calculated tangential velocity")
            sys.stdout.flush()
            
            rel_kep_tang = tang_vel/keplerian_velocity
            del tang_vel, keplerian_velocity
            gc.collect()
            
            disc_tang = np.where((rel_kep_tang>0.9)&(rel_kep_tang<1.1))[0]
            disc_full = np.where((rel_kep_full>0.9)&(rel_kep_full<1.1))[0]

            
            #get median radius of kep mass
            kep_rad_tang = np.median(radii[disc_tang])
            save_dict["Radius_tang"] = np.append(save_dict["Radius_tang"],kep_rad_tang)
            kep_rad_full = np.median(radii[disc_full])
            save_dict["Radius_full"] = np.append(save_dict["Radius_full"],kep_rad_full)
            del radii
            gc.collect()
            
            
            
            #Get keplerian mass
            kep_mass_tang = np.sum(gas_mass[disc_tang].in_units('msun'))
            save_dict["Mass_tang"] = np.append(save_dict["Mass_tang"],kep_mass_tang)
            kep_mass_full = np.sum(gas_mass[disc_full].in_units('msun'))
            save_dict["Mass_full"] = np.append(save_dict["Mass_full"],kep_mass_full)
            del gas_mass
            gc.collect
            
            #Save BHL Calculation
            file_open = open('Kep_mass_'+str(proj_root_rank)+'.pkl', 'wb')
            #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
            pickle.dump((save_dict), file_open)
            file_open.close()
            print("RANK "+str(rank)+": updated pickle", fn)
            sys.stdout.flush()

print('Finished BHL Calculation on rank', rank)
CW.Barrier()

if rank == 0:
    pickle_files = sorted(glob.glob("Kep_mass_*.pkl"))
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"Radius_tang": np.array([])})
    save_dict.update({"Radius_full": np.array([])})
    save_dict.update({"Mass_tang": np.array([])})
    save_dict.update({"Mass_full": np.array([])})
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        save_dict_r = pickle.load(file_open)
        file_open.close()
        for key in save_dict_r.keys():
            save_dict[key] = np.append(save_dict[key], save_dict_r[key])
    del save_dict_r
    gc.collect()
    sorted_inds = np.argsort(save_dict["Time"])
    for key in save_dict.keys():
        save_dict[key] = save_dict[key][sorted_inds]
        
    file_open = open('Kep_mass.pkl', 'wb')
    pickle.dump((save_dict), file_open)
    file_open.close()
    
    import matplotlib.pyplot as plt
    
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
    })

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    #Ploting parameters
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 9

    
    plt.clf()
    fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
    #axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
    plt.plot(save_dict['Time'], save_dict['Radius_tang'], label="Tangential")
    plt.plot(save_dict['Time'], save_dict['Radius_full'], label="Full")
    plt.xlabel("Time (yr)")
    plt.ylabel("$Radius (au)$")
    #plt.ylim([0, 2])
    plt.xlim([0, save_dict['Time'][-1]])
    plt.legend()
    #plt.axhline(y=0.8, ls="--", c='k')
    #plt.axhline(y=1.2, ls="--", c='k')
    #cb = fig.colorbar(smap)
    plt.savefig("Kep_mass_radius.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('Saved figure with BHL Accretion')
    
