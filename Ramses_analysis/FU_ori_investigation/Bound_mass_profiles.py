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
import matplotlib.pyplot as plt

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

parser = argparse.ArgumentParser()
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=10000, type=float)
parser.add_argument("-event_id", "--event_identifier", default=None, type=int)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------

#Start by loading pickel data and then deleting what we don't need
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

if args.event_identifier == None:
    sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/data/'
else:
    event_it = args.event_identifier
    sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/'
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
    for fn in yt.parallel_objects(files, njobs=size/28):
        frame_name = "Profile_frame" + ("%06d" % (files.index(fn)))
        if os.path.exists(frame_name+".pkl"):
            skip = True
        else:
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
            sys.stdout.flush()
            
            #Get inds in measuring sphere
            dx = ds.r['ramses', 'x'].in_units('au') - sink_pos[0].in_units('au')
            dy = ds.r['ramses', 'y'].in_units('au') - sink_pos[1].in_units('au')
            dz = ds.r['ramses', 'z'].in_units('au') - sink_pos[2].in_units('au')
            sep = np.sqrt(dx**2 + dy**2 + dz**2)
            del dx, dy, dz, sink_pos
            gc.collect()
            sys.stdout.flush()
            
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
            profile_dict = {}
            #profile_dict.update({"R_profile_mean":np.array([])})
            #profile_dict.update({"R_profile_std":np.array([])})
            #profile_dict.update({"E_profile_mean":np.array([])})
            #profile_dict.update({"E_profile_std":np.array([])})
            prev_enclosed_gas_mass = yt.YTQuantity(0, "msun")
            prev_radius = 0
            for sto, radius_bit in yt.parallel_objects(range(1, len(radius_bins)), storage=profile_dict, njobs=int(size/14)):
                #for radius_bit in range(1, len(radius_bins)):
                #Calculate enclosed mass:
                print("Calculating boundness for shell radius", radius_bins[radius_bit], "on rank", rank)
                
                enclosed_inds = np.where((sep>prev_radius)&(sep<=radius_bins[radius_bit]))[0]
                enclosed_mass = np.sum(ds.r["gas", "mass"][enclosed_inds]) + prev_enclosed_gas_mass
                prev_enclosed_gas_mass = enclosed_mass
                prev_radius = radius_bins[radius_bit]
                del enclosed_inds
                gc.collect()
                enclosed_sinks = np.where(sink_separations<=radius_bins[radius_bit])[0]
                enclosed_sink_mass = np.sum(ds.r["gas", "sink_particle_mass"][enclosed_sinks])
                del enclosed_sinks
                gc.collect()
                enclosed_mass = enclosed_mass + enclosed_sink_mass
                del enclosed_sink_mass
                gc.collect()
                
                #Now get indices in sphere
                sphere_inds = np.where((sep>radius_bins[radius_bit-1])&(sep<=radius_bins[radius_bit]))[0]
                sph_dvx = ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s') - sink_vel[0]
                sph_dvy = ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s') - sink_vel[1]
                sph_dvz = ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s') - sink_vel[2]
                sph_vel = np.sqrt(sph_dvx**2 + sph_dvy**2 + sph_dvz**2)
                del sph_dvx, sph_dvy, sph_dvz
                gc.collect()
                E_kin = 0.5 * ds.r["gas", "mass"][sphere_inds] * sph_vel**2
                del sph_vel
                gc.collect()
                #get average radius in the bin
                rad_mean = np.mean(sep[sphere_inds])
                rad_std = np.std(sep[sphere_inds])
                #calcualte gravitational potential energy
                E_grav = -1*(yt.units.gravitational_constant_cgs*enclosed_mass*ds.r["gas", "mass"][sphere_inds])/sep[sphere_inds]
                del sphere_inds
                gc.collect()
                E_ratio = E_grav.in_units('erg')/E_kin.in_units('erg')
                del E_grav, E_kin
                gc.collect()
                E_ratio_mean = np.mean(E_ratio)
                E_ratio_std = np.std(E_ratio)
                sto.result_id = str(radius_bins[radius_bit])
                sto.result = np.array([rad_mean, rad_std, E_ratio_mean, E_ratio_std])
            
            root_rank = int(rank/28)
            if rank == root_rank:
                #sort profile data:
                Profile_rad_mean = np.array([])
                Profile_rad_std = np.array([])
                Profile_E_ratio_mean = np.array([])
                Profile_E_ratio_std = np.array([])
                for key in profile_dict.keys():
                    Profile_rad_mean = np.append(Profile_rad_mean, profile_dict[key][0])
                    Profile_rad_std = np.append(Profile_rad_std, profile_dict[key][1])
                    Profile_E_ratio_mean = np.append(Profile_E_ratio_mean, profile_dict[key][2])
                    Profile_E_ratio_std = np.append(Profile_E_ratio_std, profile_dict[key][3])
                #sort inds
                sorted_inds = np.argsort(Profile_rad_mean)
                profile_dict = {}
                profile_dict.update({"R_profile_mean":Profile_rad_mean[sorted_inds]})
                profile_dict.update({"R_profile_std":Profile_rad_std[sorted_inds]})
                profile_dict.update({"E_profile_mean":Profile_E_ratio_mean[sorted_inds]})
                profile_dict.update({"E_profile_std":Profile_E_ratio_std[sorted_inds]})
            
                #Save BHL Calculation
                file_open = open(frame_name+'.pkl', 'wb')
                #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
                pickle.dump((profile_dict), file_open)
                file_open.close()
                print("RANK "+str(rank)+": updated pickle", fn)
                sys.stdout.flush()
                
                #Radial profile calcaluated, so now let's plot the frame!
                plt.clf()
                plt.xscale("log")
                plt.errorbar(profile_dict["R_profile_mean"], profile_dict["E_profile_mean"], xerr=profile_dict["R_profile_std"], yerr=profile_dict["E_profile_std"])
                plt.xlabel("Radius (au)")
                plt.ylabel("E_grav/E_kin")
                plt.xlim([np.min(profile_dict["R_profile_mean"]), np.max(profile_dict["R_profile_mean"])])
                plt.axhline(y=1.0)
                plt.savefig(frame_name+".png")

print('Finished BHL Calculation on rank', rank)
CW.Barrier()
