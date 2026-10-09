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
div_procs = 28


if len(files)>0:
    for fn in yt.parallel_objects(files, njobs=int(size/div_procs)):
        root_rank = int(rank/div_procs)
        frame_name = "Scatter_frame_" + ("%06d" % (files.index(fn)))
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
            
            #Let's doa really simple analysis ingorning enclosed gas
            sink_mass = ds.r["gas", "sink_particle_mass"][sink_id]
            
            #Get sink position
            sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
            sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
            sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
            sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
            
            #Get inds in measuring sphere
            dx = ds.r['ramses', 'x'].in_units('au') - sink_pos[0].in_units('au')
            dy = ds.r['ramses', 'y'].in_units('au') - sink_pos[1].in_units('au')
            dz = ds.r['ramses', 'z'].in_units('au') - sink_pos[2].in_units('au')
            sep = np.sqrt(dx**2 + dy**2 + dz**2)
            del dx, dy, dz, sink_pos
            gc.collect()
            sys.stdout.flush()
            
            usuable_inds = np.where(sep<10000)[0]
            sep = sep[usuable_inds]

            E_grav = -1*(yt.units.gravitational_constant_cgs*sink_mass*ds.r["gas", "mass"][usuable_inds])/sep
            
            
            sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
            sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
            sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
            sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
            del sink_particle_velx, sink_particle_vely, sink_particle_velz
            gc.collect()
            print("RANK", rank, "got particle velocity")
            
            dvx = ds.r["ramses", "x-velocity"][usuable_inds].in_units('km/s') - sink_vel[0]
            dvy = ds.r["ramses", "y-velocity"][usuable_inds].in_units('km/s') - sink_vel[1]
            dvz = ds.r["ramses", "z-velocity"][usuable_inds].in_units('km/s') - sink_vel[2]
            vel = np.sqrt(dvx**2 + dvx**2 + dvx**2)
            del dvx, dvy, dvz, sink_vel
            gc.collect()
            sys.stdout.flush()
            
            E_kin = 0.5 * ds.r["gas", "mass"][usuable_inds] * vel**2
            
            E_tot = E_grav.in_units('erg')+E_kin.in_units('erg')
            bound_inds = np.argwhere(E_tot<0)
            
            if rank == root_rank:
                #Radial profile calcaluated, so now let's plot the frame!
                plt.clf()
                plt.scatter(sep, abs(E_kin.in_units('erg'))/abs(E_grav.in_units('erg')), c=ds.r["gas", "mass"][usuable_inds], marker='.', rasterized=True)
                plt.xlabel('Separation (AU)')
                plt.ylabel('E_kin/E_grav')
                plt.yscale("log")
                plt.axhline(y=1.0)
                plt.savefig(frame_name+".png")
                print("saved scatter plot on rank,")
                plt.clf()
                gc.collect()
            

print('Finished BHL Calculation on rank', rank)
CW.Barrier()
