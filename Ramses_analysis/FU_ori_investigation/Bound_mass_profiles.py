#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import os
import yt
import my_ramses_fields_short as myf
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
parser.add_argument("-para_div", "--parallelisation_jobs", default=None, type=int)
parser.add_argument("-use_com", "--use_centre_of_mass", default=False, type=str)
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

sink_id = 45
sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()
gc.collect()

#Define radius bins
radius = yt.YTQuantity(args.measuring_sphere_radius, 'au')
radius_bins = np.logspace(0, np.log10(radius), 50)
if args.parallelisation_jobs == None:
    para_div = size
else:
    para_div = args.parallelisation_jobs

if len(files)>0:
    for fn in yt.parallel_objects(files, njobs=size):
        frame_name = "Profile_frame_" + ("%06d" % (files.index(fn)))
        if os.path.exists(frame_name+".pkl"):
            skip = True
        else:
            print('Reading file', fn, 'on rank', rank)
            sys.stdout.flush()
            ds = yt.load(fn, units_override=units_override)
            sinks = myf.get_sink_table(ds)
            
            if len(sinks['mass']) == 45:
                skip=True
            else:
                if np.isnan(sink_form_time):
                    sink_form_time = sinks['form_time'][sink_id]
                    print("RANK", rank, "got sink formation time")
                    sys.stdout.flush()
                skip = False
        if skip == False:
            time_val = ds.current_time.in_units('yr').value - sink_form_time.in_units('yr').value
            gc.collect()
            sys.stdout.flush()
            #Get sink position
            
            sinks = myf.get_sink_table(ds)
            center_pos = sinks['pos'][sink_id]
            sink_separations = np.sqrt(np.sum((sinks['pos'] - center_pos)**2, axis=1))
            
            if args.use_centre_of_mass != False:
                nearest_sink = np.argsort(sink_separations)[1]
                binary_masses = sinks['mass'][[sink_id, nearest_sink]]
                binary_positions = sinks['pos'][[sink_id, nearest_sink]]
                binary_velocities = sinks['vel'][[sink_id, nearest_sink]]
                center_pos = np.sum(binary_masses*binary_positions.T, axis=1)/np.sum(binary_masses)
                center_vel = np.sum(binary_masses*binary_velocities.T, axis=1)/np.sum(binary_masses)
                sink_separations = np.sqrt(np.sum((sinks['pos'] - center_pos)**2, axis=1))
            else:
                center_vel = sinks['vel'][sink_id]
            
            sink_masses = sinks['mass']
            sink_masses = np.delete(sink_masses, np.where(sink_separations>radius)[0])
            sink_separations = np.delete(sink_separations, np.where(sink_separations>radius)[0])
            #Delete sink separations greater than measuring sphere
            gc.collect()
            sys.stdout.flush()
            
            #Get inds in measuring sphere
            dx = ds.r['ramses', 'x'].in_units('au') - center_pos[0].in_units('au')
            dy = ds.r['ramses', 'y'].in_units('au') - center_pos[1].in_units('au')
            dz = ds.r['ramses', 'z'].in_units('au') - center_pos[2].in_units('au')
            sep = np.sqrt(dx**2 + dy**2 + dz**2)
            del dx, dy, dz, center_pos
            del_inds = np.where(sep>radius)[0]
            gas_mass = ds.r["gas", "mass"]
            gas_mass = np.delete(gas_mass, np.where(sep>radius)[0])
            sep = np.delete(sep, np.where(sep>radius)[0])
            #delete cell further than measuring sphere
            gc.collect()
            sys.stdout.flush()
            
            
            gc.collect()
            print("RANK", rank, "got particle velocity")
            sys.stdout.flush()
            
            sph_dvx = ds.r["ramses", "x-velocity"].in_units('km/s') - center_vel[0]
            sph_dvy = ds.r["ramses", "y-velocity"].in_units('km/s') - center_vel[1]
            sph_dvz = ds.r["ramses", "z-velocity"].in_units('km/s') - center_vel[2]
            sph_vel = np.sqrt(sph_dvx**2 + sph_dvy**2 + sph_dvz**2)
            sph_vel = np.delete(sph_vel, del_inds)
            del del_inds
            
            #get separations to other sink particles
            
            #Get indices in measure sphere
            #Start iterating over Radial bins
            profile_dict = {}
            profile_dict.update({"time":time_val})
            profile_dict.update({"R_profile_mean":np.array([])})
            profile_dict.update({"R_profile_std":np.array([])})
            profile_dict.update({"E_ratio_profile_mean":np.array([])})
            profile_dict.update({"E_ratio_profile_std":np.array([])})
            profile_dict.update({"E_tot_profile_mean":np.array([])})
            profile_dict.update({"E_tot_profile_std":np.array([])})
            profile_dict.update({"Sink_seps":np.array([])})
            prev_enclosed_gas_mass = yt.YTQuantity(0, "msun")
            prev_radius = 0
            #for sto, radius_bit in yt.parallel_objects(range(1, len(radius_bins)), storage=profile_dict, njobs=1):
            for radius_bit in range(1, len(radius_bins)):
                #for radius_bit in range(1, len(radius_bins)):
                #Calculate enclosed mass:
                
                del_inds = np.where(sep<prev_radius)[0]
                sep = np.delete(sep, del_inds)
                sph_vel = np.delete(sph_vel, del_inds)
                gc.collect()
                enclosed_inds = np.where((sep>prev_radius)&(sep<=radius_bins[radius_bit]))[0]
                enclosed_mass = np.sum(gas_mass[enclosed_inds]) + prev_enclosed_gas_mass
                #FIGURE OUT WHY ENCLOSED MASS IC GROWING CRAZY FAST
                prev_enclosed_gas_mass = enclosed_mass
                prev_radius = radius_bins[radius_bit]
                del enclosed_inds
                gc.collect()
                enclosed_sinks = np.where(sink_separations<=radius_bins[radius_bit])[0]
                enclosed_sink_mass = np.sum(sink_masses[enclosed_sinks])
                del enclosed_sinks
                gc.collect()
                enclosed_mass = enclosed_mass + enclosed_sink_mass
                del enclosed_sink_mass
                gc.collect()
                
                #Now get indices in sphere
                sphere_inds = np.where((sep>radius_bins[radius_bit-1])&(sep<=radius_bins[radius_bit]))[0]
                gc.collect()
                E_kin = 0.5 * gas_mass[sphere_inds] * sph_vel[sphere_inds]**2
                gc.collect()
                #get average radius in the bin
                rad_mean = np.mean(sep[sphere_inds])
                profile_dict["R_profile_mean"] = np.append(profile_dict["R_profile_mean"], rad_mean)
                rad_std = np.std(sep[sphere_inds])
                profile_dict["R_profile_std"] = np.append(profile_dict["R_profile_std"], rad_std)
                #calcualte gravitational potential energy
                E_grav = -1*(yt.units.gravitational_constant_cgs*enclosed_mass*gas_mass[sphere_inds])/sep[sphere_inds]
                del sphere_inds
                gc.collect()
                E_tot = E_grav.in_units('erg')+E_kin.in_units('erg')
                E_ratio = (-1*E_grav.in_units('erg'))/E_kin.in_units('erg')
                del E_grav, E_kin
                gc.collect()
                E_tot_mean = np.mean(E_tot)
                E_tot_std = np.std(E_tot)
                E_ratio_mean = np.mean(E_ratio)
                E_ratio_std = np.std(E_ratio)
                
                profile_dict["E_ratio_profile_mean"] = np.append(profile_dict["E_ratio_profile_mean"], E_ratio_mean)
                profile_dict["E_ratio_profile_std"] = np.append(profile_dict["E_ratio_profile_std"], E_ratio_std)
                profile_dict["E_tot_profile_mean"] = np.append(profile_dict["E_tot_profile_mean"], E_tot_mean)
                profile_dict["E_tot_profile_std"] = np.append(profile_dict["E_tot_profile_std"], E_tot_std)
                profile_dict["Sink_seps"] = np.append(profile_dict["Sink_seps"], sink_separations)
                del E_ratio, E_tot
                gc.collect()
                #sto.result_id = str(radius_bins[radius_bit])
                #sto.result = np.array([rad_mean, rad_std, E_ratio_mean, E_ratio_std, E_tot_mean, E_tot_std])
                print("Calculated boundness for shell radius", radius_bins[radius_bit], "on rank", rank, ": Enclosed_mass, E_ratio_mean, E_tot_mean", enclosed_mass.in_units('msun'), E_ratio_mean, E_tot_mean)
                sys.stdout.flush()
            
            #Save BHL Calculation
            file_open = open(frame_name+'.pkl', 'wb')
            #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
            pickle.dump((profile_dict), file_open)
            file_open.close()
            gc.collect()
            print("RANK "+str(rank)+": updated pickle", fn)
            sys.stdout.flush()
                
        if os.path.exists(frame_name+".pkl") and rank == 0:
            file_open = open(frame_name+'.pkl', 'rb')
            #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
            profile_dict = pickle.load(file_open)
            file_open.close()
                
            #Radial profile calcaluated, so now let's plot the frame!
            plt.clf()
            plt.xscale("log")
            plt.title("time =" + str(profile_dict["time"]))
            plt.errorbar(profile_dict["R_profile_mean"], profile_dict["E_ratio_profile_mean"], xerr=profile_dict["R_profile_std"], yerr=profile_dict["E_ratio_profile_std"])
            for sink_sep in profile_dict["Sink_seps"]:
                plt.axvline(x=sink_sep, c='r')
            plt.axhline(y=1.0, ls='--', c='k')
            plt.xlabel("Radius (au)")
            plt.ylabel("E_grav/E_kin")
            plt.xlim([np.min(profile_dict["R_profile_mean"]), np.max(profile_dict["R_profile_mean"])])
            plt.ylim([0.0, np.max(profile_dict["E_ratio_profile_mean"])])
            plt.savefig(frame_name+"_ratio.jpg")
            
            plt.clf()
            plt.xscale("log")
            plt.title("time =" + str(profile_dict["time"]))
            plt.errorbar(profile_dict["R_profile_mean"], profile_dict["E_tot_profile_mean"], xerr=profile_dict["R_profile_std"], yerr=profile_dict["E_tot_profile_std"])
            for sink_sep in profile_dict["Sink_seps"]:
                plt.axvline(x=sink_sep)
            plt.xlabel("Radius (au)")
            plt.ylabel("E_grav/E_kin")
            plt.xlim([np.min(profile_dict["R_profile_mean"]), np.max(profile_dict["R_profile_mean"])])
            plt.ylim([np.min(profile_dict["E_tot_profile_mean"]), np.max(profile_dict["E_tot_profile_mean"])])
            plt.yscale("symlog")
            plt.savefig(frame_name+"_total.jpg")
            del profile_dict
            
            plt.clf()
        gc.collect()
        CW.Barrier()

print('Finished BHL Calculation on rank', rank)
CW.Barrier()
