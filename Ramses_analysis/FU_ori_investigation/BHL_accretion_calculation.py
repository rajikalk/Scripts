#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import matplotlib.pyplot as plt
import os
import yt
yt.enable_parallelism()
import my_ramses_fields_short as myf
import gc
from mpi4py.MPI import COMM_WORLD as CW

def projected_vector(vector, proj_vector):
    """
    Calculates the position of vector projected onto proj_vector
    """
    vector_units = vector.units
    proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
    proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
    proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
    proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-ax", "--axis", default='xy', type=str)
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=6, type=float)
parser.add_argument('files', nargs='*')
args = parser.parse_args()


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

#------------------------------------------------------
time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]
cmap=plt.cm.gist_heat

#Start by loading pickel data and then deleting what we don't need

event_it = args.event_identifier
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

sys.stdout.flush()
CW.Barrier()

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/'
files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))#[::10]


if os.path.exists('BHL_accretion.pkl'):
    file_open = open('BHL_accretion.pkl', 'rb')
    save_dict = pickle.load(file_open)
    file_open.close()
    if len(save_dict["Time"]) != len(files):
        files = files[len(save_dict["Time"]):]
elif os.path.exists('BHL_accretion_0.pkl'):
    pickle_files = sorted(glob.glob("BHL_accretion_*.pkl"))
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"BHL_Acc_acc_low": np.array([])})
    save_dict.update({"BHL_Acc_acc_high": np.array([])})
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
        files = files[len(time_arr):]
else:
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"BHL_Acc_acc_low": np.array([])})
    save_dict.update({"BHL_Acc_acc_high": np.array([])})
    save_dict.update({"Density": np.array([])})
    save_dict.update({"Rel_kep": np.array([])})

sink_id = 45
sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()

#sink_form_time = np.nan

sys.stdout.flush()
CW.Barrier()

if len(files)>0:
    #ts = yt.DatasetSeries(files, parallel=4)
    #'''
    para_div = 7
    #my_storage = {}
    for fn in yt.parallel_objects(files, njobs=int(size/para_div)):#, storage=my_storage):
        proj_root_rank = int(rank/(size/para_div))
        print('Reading file', fn, 'on rank', rank)
        sys.stdout.flush()
        #try:
        ds = yt.load(fn, units_override=units_override)
        #'''
        #my_storage = {}
        #for sto, ds in ts.piter(storage=my_storage):
        if np.isnan(sink_form_time):
            sink_form_time = ds.r["sink_particle_form_time"][sink_id]
        time_val = ds.current_time.in_units('yr').value - sink_form_time.in_units('yr').value
        save_dict["Time"] = np.append(save_dict["Time"],time_val)
        #sto.result_id = "Time"
        #sto.result = time_val
        
        sink_mass = ds.r["gas", "sink_particle_mass"][sink_id]
        gc.collect()
        print('Got particle mass on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        #Get sink position
        sink_particle_posx = ds.r["gas", "sink_particle_posx"][sink_id]
        sink_particle_posy = ds.r["gas", "sink_particle_posy"][sink_id]
        sink_particle_posz = ds.r["gas", "sink_particle_posz"][sink_id]
        sink_pos = yt.YTArray([sink_particle_posx, sink_particle_posy, sink_particle_posz])
        del sink_particle_posx, sink_particle_posy, sink_particle_posz
        gc.collect()
        print('Got particle position on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        #get sink velocity
        sink_particle_velx = ds.r["gas", "sink_particle_velx"][sink_id]
        sink_particle_vely = ds.r["gas", "sink_particle_vely"][sink_id]
        sink_particle_velz = ds.r["gas", "sink_particle_velz"][sink_id]
        sink_vel = yt.YTArray([sink_particle_velx, sink_particle_vely, sink_particle_velz])
        del sink_particle_velx, sink_particle_vely, sink_particle_velz
        gc.collect()
        print('Got particle velocity on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        #Define measuring sphere:
        radius = yt.YTQuantity(args.measuring_sphere_radius, 'au')
        
        #Get inds in measuring sphere
        dd = ds.all_data()
        dx = dd['x'].in_units('au') - sink_pos[0].in_units('au')
        dy = dd['y'].in_units('au') - sink_pos[1].in_units('au')
        dz = dd['z'].in_units('au') - sink_pos[2].in_units('au')
        sep_vector = yt.YTArray([dx, dy, dz])
        del dx, dy, dz, dd
        gc.collect()
        sep = np.sqrt(sep_vector[0]**2 + sep_vector[1]**2 + sep_vector[2]**2)
        print('Got indexes of cells in measuring sphere on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        #Get indices in measure sphere
        sphere_inds = np.where(sep<=radius)[0]
        sep_vector = sep_vector.T[sphere_inds].T
        
        #Calcualte keplerian velocity
        radii = sep[sphere_inds]
        del sep
        gc.collect()
        gas_mass = ds.r["gas", "mass"][sphere_inds]
        enclosed_mass = yt.YTArray(np.zeros(np.shape(radii)), "g")
        for radi_it in range(len(radii)):
            enc_inds = np.where(radii<=radii[radi_it])[0]
            enc_mass = np.sum(gas_mass[enc_inds])
            enclosed_mass[radi_it] = enc_mass
        del gas_mass
        gc.collect()
        enclosed_mass = enclosed_mass+sink_mass.in_units('g')
        keplerian_velocity = np.sqrt((yt.units.gravitational_constant_cgs*enclosed_mass)/radii).in_units('km/s')
        
        density_array = ds.r["gas", "Density"][sphere_inds]
        save_dict["Density"] = np.append(save_dict["Density"], density_array)
        mean_density = np.mean(density_array)
        del density_array
        gc.collect()
        
        #Calculate bulk velocity of the sphere
        sph_dvx = ds.r["ramses", "x-velocity"][sphere_inds].in_units('km/s') - sink_vel[0]
        sph_dvy = ds.r["ramses", "y-velocity"][sphere_inds].in_units('km/s') - sink_vel[1]
        sph_dvz = ds.r["ramses", "z-velocity"][sphere_inds].in_units('km/s') - sink_vel[2]
        sph_vel = yt.YTArray([sph_dvx, sph_dvy, sph_dvz])
        rel_vel = yt.YTArray([np.mean(sph_dvx), np.mean(sph_dvy), np.mean(sph_dvz)])
        del sph_dvx, sph_dvy, sph_dvz
        gc.collect()
        
        #Calculate Tangential vel
        proj_v_x = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[0]
        proj_v_y = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[1]
        proj_v_z = (np.dot(sph_vel.T.in_units('km/s'), sep_vector.in_units('km')).diagonal())/np.dot(sep_vector.T.in_units('km'), sep_vector.in_units('km')).diagonal()*sep_vector.in_units('km')[2]
        del sep_vector
        gc.collect()
        rad_vel = yt.YTArray([proj_v_x,proj_v_y,proj_v_z])
        del proj_v_x, proj_v_y, proj_v_z
        gc.collect()
        rad_speed = np.sqrt(rad_vel[0]**2 + rad_vel[1]**2 + rad_vel[2]**2)
        del rad_vel
        gc.collect()
        sph_speed = np.sqrt(sph_vel[0]**2 + sph_vel[1]**2 + sph_vel[2]**2)
        del sph_vel
        gc.collect()
        tang_vel = np.sqrt(sph_speed**2 - rad_speed**2)
        del sph_speed, rad_speed
        gc.collect()
        rel_kep = tang_vel/keplerian_velocity
        save_dict["Rel_kep"] = np.append(save_dict["Rel_kep"], rel_kep)
        del rel_kep
        gc.collect()
        
        rel_speed = np.sqrt(np.sum(rel_vel**2))
        del rel_vel
        gc.collect()
        print('calculated mean density and relative speed on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        sound_speed = np.mean(np.sqrt((ds.r["gas", "Gamma"][sphere_inds]*ds.r["gas", "Pressure"][sphere_inds])/ds.r["gas", "Density"][sphere_inds]).in_units('km/s'))
        print('calculated sound speed on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        alpha = yt.YTArray([1, 2], '')
        BHL_top = 2*np.pi* (sink_mass.in_cgs()*yt.units.gravitational_constant_cgs)**2 * mean_density.in_cgs()
        del sink_mass, mean_density
        gc.collect()
        BHL_bot = (rel_speed.in_cgs()**2 + sound_speed.in_cgs()**2)**(3./2.)
        del rel_speed, sound_speed
        gc.collect()
        BHL = (alpha * (BHL_top/BHL_bot)).in_units('msun/yr')
        del alpha, BHL_top, BHL_bot
        gc.collect()
        save_dict["BHL_Acc_acc_low"] = np.append(save_dict["BHL_Acc_acc_low"], BHL[0])
        save_dict["BHL_Acc_acc_high"] = np.append(save_dict["BHL_Acc_acc_high"], BHL[1])
        #sto.result_id = "BHL_Acc_acc_low"
        #sto.result = BHL[0]
        #sto.result_id = "BHL_Acc_acc_high"
        #sto.result = BHL[1]
        print('calculated BHL accretion on rank', rank, ' for fn', ds)
        sys.stdout.flush()
        
        #ds.index.clear_all_data()
        
        #Save BHL Calculation
        file_open = open('BHL_accretion_'+str(proj_root_rank)+'.pkl', 'wb')
        #pickle.dump((my_storage["Time"], my_storage["BHL_Acc_acc_low"], my_storage["BHL_Acc_acc_high"]), file_open)
        pickle.dump((save_dict), file_open)
        file_open.close()
        print("RANK "+str(rank)+": Calculated BHL for file", fn)
        sys.stdout.flush()
            
        #except:
        #    print(fn, "seems to be missing data")


print('Finished BHL Calculation on rank', rank)
CW.Barrier()

if rank == 0:
    pickle_files = sorted(glob.glob("BHL_accretion_*.pkl"))
    save_dict = {}
    save_dict.update({"Time": np.array([])})
    save_dict.update({"BHL_Acc_acc_low": np.array([])})
    save_dict.update({"BHL_Acc_acc_high": np.array([])})
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
    
    try:
        sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
        file_open = open(sink_pickle, 'rb')
        particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
        file_open.close()
        print("finished reading in pickle")
        sys.stdout.flush()
    except:
        sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
        print("read pickle", sink_pickle)
        file_open = open(sink_pickle, 'rb')
        particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
        file_open.close()
        print("finished reading in pickle")
        sys.stdout.flush()
    
    plt.clf()
    fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
                
    plt.title("Burst event "+str(event_it), y=0.8)
    start_ind = np.argmin(abs(particle_data['time']-start_time))
    end_ind = np.argmin(abs(particle_data['time']-end_time))
    #axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
    lns1 = plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='b', ls='-', label="Accretion rate")
    BHL_mean = (save_dict["BHL_Acc_acc_low"]+save_dict["BHL_Acc_acc_high"])/2
    lns3 = plt.semilogy(save_dict["Time"], BHL_mean, color='g', ls=':', label="BHL_mean")
    plt.fill_between(time_arr, save_dict["BHL_Acc_acc_low"], save_dict["BHL_Acc_acc_high"], color='g', alpha=0.5, label="BHL prediction")
    axes_1_twin = plt.twinx()
    lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
    #Plot accretion and separation. This should be loaded from a pickle

    plt.xlabel('Time (yr)', labelpad=-0.2, fontsize=font_size) #($yr$)
    plt.ylabel('Accretion rate (M$_\odot$/yr)', labelpad=-0.2, fontsize=font_size)# (M$_\odot/yr$)
    axes_1_twin.set_ylabel('Separation (au)', fontsize=font_size)
    plt.tick_params(axis='x', which='major', direction='in', color='k', top=True)
    plt.tick_params(axis='y', which='major', direction='in', color='k', right=True)
    #plt.xaxis.label.set_color('black')
    #plt.yaxis.label.set_color('black')
    plt.tick_params(axis='both', labelsize=font_size)
    plt.xlim([start_time, end_time])
    plt.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='upper left')
    plt.savefig("BHL_Event_"+str(event_it)+"_radius_"+str(args.measuring_sphere_radius)+"au.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('Saved figure with BHL Accretion')

    #plt.savefig("BHL_Event_"+str(event_it)+".pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='upper left')
    plt.savefig("BHL_Event_"+str(event_it)+"_radius_"+str(args.measuring_sphere_radius)+"au.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('Saved figure with BHL Accretion')
    sys.stdout.flush()
    
