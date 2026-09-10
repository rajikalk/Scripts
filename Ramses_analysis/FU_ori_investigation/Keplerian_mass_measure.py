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
import my_ramses_fields_short as myf
import gc
from mpi4py.MPI import COMM_WORLD as CW

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

parser = argparse.ArgumentParser()
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=10, type=float)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------

#Start by loading pickel data and then deleting what we don't need
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

sys.stdout.flush()
CW.Barrier()

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/data/'
files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))#[::10]

if os.path.exists('Kep_mass.pkl'):
    file_open = open('Kep_mass.pkl', 'rb')
    save_dict = pickle.load(file_open)
    file_open.close()
    if len(save_dict["Time"]) != len(files):
        files = files[len(save_dict["Time"]):]
elif os.path.exists('Kep_mass.pkl'):
    pickle_files = sorted(glob.glob("Kep_mass_*.pkl"))
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"Density": np.array([])})
    save_dict.update({"Rel_kep": np.array([])})
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        save_dict_r = pickle.load(file_open)
        file_open.close()
        for key in save_dict_r.keys():
            save_dict[key] = np.append(save_dict[key],save_dict_r[key])
    sorted_inds = np.argsort(save_dict["Time"])
    for key in save_dict.keys():
        save_dict[key] = save_dict[key][sorted_inds]
    if len(save_dict["Time"]) != len(files):
        files = files[len(save_dict["Time"]):]
else:
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"Density": np.array([])})
    save_dict.update({"Rel_kep": np.array([])})

sink_id = 45
sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()
gc.collect()

if len(files)>0:
    #ts = yt.DatasetSeries(files, parallel=4)
    #'''
    para_div = 7
    #my_storage = {}
    for fn in yt.parallel_objects(files, njobs=int(size/para_div)):#, storage=my_storage):
        proj_root_rank = int(rank/para_div)
        print('Reading file', fn, 'on rank', rank)
        sys.stdout.flush()
        #try:
        ds = yt.load(fn, units_override=units_override)
        #'''
        #my_storage = {}
        #for sto, ds in ts.piter(storage=my_storage):
        if np.isnan(sink_form_time):
            if len(ds.r["sink_particle_form_time"]) == 45:
                skip=True
            else:
                sink_form_time = ds.r["sink_particle_form_time"][sink_id]
                skip = False
        if skip == False:
            time_val = ds.current_time.in_units('yr').value - sink_form_time.in_units('yr').value
            save_dict["Time"] = np.append(save_dict["Time"],time_val)
            del time_val
            gc.collect()
            #sto.result_id = "Time"
            #sto.result = time_val
            
            #Get sink position
            sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
            sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
            sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
            sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
            del sink_particle_posx, sink_particle_posy, sink_particle_posz
            gc.collect()
            #print('Got particle position on rank', rank, ' for fn', ds)
            sys.stdout.flush()
            
            #Get inds in measuring sphere
            dd = ds.all_data()
            dx = dd['x'].in_units('au') - sink_pos[0].in_units('au')
            dy = dd['y'].in_units('au') - sink_pos[1].in_units('au')
            dz = dd['z'].in_units('au') - sink_pos[2].in_units('au')
            sep_vector = yt.YTArray([dx, dy, dz])
            del dx, dy, dz, dd, sink_pos
            gc.collect()
            sep = np.sqrt(sep_vector[0]**2 + sep_vector[1]**2 + sep_vector[2]**2)
            sys.stdout.flush()
            
            #Get indices in measure sphere
            radius = yt.YTQuantity(args.measuring_sphere_radius, 'au')
            sphere_inds = np.where(sep<=radius)[0]
            radii = sep[sphere_inds]
            del sep
            gc.collect()
            #print('Got indexes of cells in measuring sphere on rank', rank, ' for fn', ds)
            sep_vector = sep_vector.T[sphere_inds].T
            
            #Calculate bulk velocity of the sphere
            sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
            sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
            sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
            sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
            del sink_particle_velx, sink_particle_vely, sink_particle_velz
            gc.collect()
            #print('Got particle velocity on rank', rank, ' for fn', ds)
            sys.stdout.flush()
            sph_dvx = ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s') - sink_vel[0]
            sph_dvy = ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s') - sink_vel[1]
            sph_dvz = ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s') - sink_vel[2]
            sph_vel = yt.YTArray([sph_dvx, sph_dvy, sph_dvz])
            del sph_dvx, sph_dvy, sph_dvz, sink_vel
            gc.collect()
            sph_speed = np.sqrt(sph_vel[0]**2 + sph_vel[1]**2 + sph_vel[2]**2)
            rel_vel = yt.YTArray([np.mean(sph_vel[0]), np.mean(sph_vel[1]), np.mean(sph_vel[2])])
            del rel_vel
            gc.collect()
            print('calculated mean density and relative speed on rank', rank, ' for fn', ds)
            sys.stdout.flush()
            
            #Calculate Tangential vel
            proj_v_x = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[0]
            proj_v_y = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[1]
            proj_v_z = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[2]
            del sph_vel, sep_vector
            gc.collect()
            rad_vel = yt.YTArray([proj_v_x,proj_v_y,proj_v_z])
            del proj_v_x, proj_v_y, proj_v_z
            gc.collect()
            rad_speed = np.sqrt(rad_vel[0]**2 + rad_vel[1]**2 + rad_vel[2]**2)
            del rad_vel
            gc.collect()
            tang_vel = np.sqrt(sph_speed**2 - rad_speed**2)
            del sph_speed, rad_speed
            gc.collect()
            
            #Calcualte keplerian velocity
            gas_mass = ds.r["gas", "mass"][sphere_inds]
            enclosed_mass = yt.YTArray(np.zeros(np.shape(radii)), "g")
            for radi_it in range(len(radii)):
                enc_inds = np.where(radii<=radii[radi_it])[0]
                enc_mass = np.sum(gas_mass[enc_inds])
                enclosed_mass[radi_it] = enc_mass
            del gas_mass
            gc.collect()
            sink_mass = ds.r["gas", "sink_particle_mass"][sink_id]
            #print('Got particle mass on rank', rank, ' for fn', ds)
            sys.stdout.flush()
            enclosed_mass = enclosed_mass+sink_mass.in_units('g')
            keplerian_velocity = np.sqrt((yt.units.gravitational_constant_cgs*enclosed_mass)/radii).in_units('km/s')
            del enclosed_mass, sink_mass
            gc.collect()
            
            rel_kep = tang_vel/keplerian_velocity
            save_dict["Rel_kep"] = np.append(save_dict["Rel_kep"], rel_kep)
            del rel_kep
            gc.collect()
            
            density_array = ds.r["gas", "Density"][sphere_inds]
            save_dict["Density"] = np.append(save_dict["Density"], density_array)
            del density_array
            gc.collect()
            
            #Save BHL Calculation
            file_open = open('Kep_mass_'+str(proj_root_rank)+'.pkl', 'wb')
            #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
            pickle.dump((save_dict), file_open)
            file_open.close()
            print("RANK "+str(rank)+": Calculated Keplerian mass for file", fn)
            sys.stdout.flush()

print('Finished BHL Calculation on rank', rank)
CW.Barrier()

if rank == 0:
    pickle_files = sorted(glob.glob("BHL_accretion_*.pkl"))
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"Density": np.array([])})
    save_dict.update({"Rel_kep": np.array([])})
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        save_dict_r = pickle.load(file_open)
        file_open.close()
        for key in save_dict_r.keys():
            save_dict[key] = np.append(save_dict[key],save_dict_r[key])
    sorted_inds = np.argsort(save_dict["Time"])
    for key in save_dict.keys():
        save_dict[key] = save_dict[key][sorted_inds]
        
    file_open = open('BHL_accretion.pkl', 'wb')
    pickle.dump((save_dict), file_open)
    file_open.close()
    
    import matplotlib.pyplot as plt
    
    cmap=plt.cm.gist_heat
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
    mym.set_global_font_size(font_size)

    
    plt.clf()
    fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
    #axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
    plt.scatter(particle_data['time'], particle_data['Rel_kep'], color=particle_data['Density'])
    plt.ylim([0, 2])
    plt.axhline(y=0.8, ls="--", c='k')
    plt.axhline(y=1.2, ls="--", c='k')
    plt.savefig("Kep_mass_radius_"+str(args.measuring_sphere_radius)+"au.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('Saved figure with BHL Accretion')
    sys.stdout.flush()
    
