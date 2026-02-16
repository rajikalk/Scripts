#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-make_pickles", "--make_pickle_files", help="do you want to update pickles?", default='True', type=str)
    parser.add_argument("-radius", "--analysis_radius", type=float, default=100.)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

#Set some plot variables independant on data files
sys.stdout.flush()
CW.Barrier()

#File files
files = sorted(glob.glob(input_dir+"*/info*.txt"))

sys.stdout.flush()
CW.Barrier()

#Define units to override:
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
mym.set_units(units_override)


#find sink particle to center on and formation time
ds = yt.load(files[-1], units_override=units_override)
#try:
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_com_pos_use_gas(False)
myf.set_com_pos_use_part(True)
myf.set_com_vel_use_gas(False)
myf.set_com_vel_use_part(True)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
#find start file
usable_files = mym.find_files([0], files, sink_form_time, sink_id, verbatim=False)
start_ind = files.index(usable_files[0])
files = files[start_ind+1:]
vel_bins = np.linspace(0, 50, 26).tolist()
vel_bins.append(100)
del dd

if args.make_pickle_files == 'True':

    sink_dict = {'time':[], 'mass':[], 'mdot':[], 'max_outflow_speed':[], 'mean_density':[], 'inflow_mass':[], 'outflow_distribution':[]}

    for fn in yt.parallel_objects(files):
        pickle_file = 'burst_analysis_sink_'+str(sink_id)+'_'+str(rank)+'.pkl'
        ds = yt.load(fn, units_override=units_override)
        dd = ds.all_data()
        
        #save time, mass, and accretion
        sink_dict['time'].append(ds.current_time.in_units('yr') - sink_form_time)
        sink_dict['mass'].append(dd['sink_particle_mass'][sink_id].in_units('msun'))
        sink_dict['mdot'].append(dd['sink_particle_accretion_rate'][sink_id].in_units('msun/yr'))
        #Define box:
        #center_pos = yt.YTArray([dd['sink_particle_posx'][sink_id].in_units('au'), dd['sink_particle_posy'][sink_id].in_units('au'), dd['sink_particle_posz'][sink_id].in_units('au')])
        #center_vel = yt.YTArray([dd['sink_particle_velx'][sink_id].in_units('cm/s'), dd['sink_particle_vely'][sink_id].in_units('cm/s'), dd['sink_particle_velz'][sink_id].in_units('cm/s')])
        
        #get position and velocity of centred sink particle
        center_pos = dd['Center_Position'].in_units('au')
        center_vel = dd['Center_Velocity'].in_units('cm/s')
        
        #Use Sphere to get Angular momentum vector
        sph = ds.sphere(center_pos, (args.analysis_radius, "au"))
        
        if args.analysis_radius > 50:
            L_vec = yt.YTArray([np.sum(sph['Angular_Momentum_x']), np.sum(sph['Angular_Momentum_y']), np.sum(sph['Angular_Momentum_z'])])
            L_mag = np.sqrt(np.sum(L_vec**2))
            L_norm = L_vec/L_mag
            
            #Define cyclinder based on momentum vector
            disk = ds.disk(center_pos, L_norm, (args.analysis_radius, "au"), (args.analysis_radius, "au"))
        else:
            disk = sph
        
        #Work our which cells have velocities towards to away from the sink
        sep_vector = yt.YTArray([disk['x'].in_units('cm')-center_pos[0].in_units('cm'), disk['y'].in_units('cm')-center_pos[1].in_units('cm'), disk['z'].in_units('cm')-center_pos[2].in_units('cm')]).T
        sep_vector_length = np.sqrt(np.sum(sep_vector**2, axis=1))
        sep_vector_norm = (sep_vector.T/sep_vector_length).T
        
        gas_vel = yt.YTArray([disk['Corrected_velx'], disk['Corrected_vely'], disk['Corrected_velz']]).in_units('cm/s').T
        #gas_vel = yt.YTArray([disk['x-velocity']-center_vel[0], disk['y-velocity']-center_vel[1], disk['z-velocity']-center_vel[2]]).in_units('cm/s').T
        gas_vel_length = np.sqrt(np.sum(gas_vel**2, axis=1))
        gas_vel_norm = (gas_vel.T/gas_vel_length).T
        
        vel_dot = gas_vel_norm.T[0]*sep_vector_norm.T[0] + gas_vel_norm.T[1]*sep_vector_norm.T[1] + gas_vel_norm.T[2]*sep_vector_norm.T[2]
        
        #Quantify mass flux in disc. I guess just summing the mass of inflowing material
        inflow_inds = np.where(vel_dot<0)[0]
        inflow_mass = np.sum(disk['cell_mass'][inflow_inds].in_units('msun'))
        sink_dict['inflow_mass'].append(inflow_mass)
        sink_dict['mean_density'].append(np.mean(disk['density'][inflow_inds]))
        
        #Get a distribution of outflow velocities
        outflow_inds = np.where(vel_dot>0)[0]
        outflow_velocities = disk['Corrected_vel_mag'][outflow_inds].in_units('km/s')

        #Make histogram
        vel_hist, vel_bins = np.histogram(outflow_velocities, bins=vel_bins)
        vel_hist_norm = (vel_hist)/np.sum(vel_hist)
        sink_dict['outflow_distribution'].append(vel_hist_norm)
        
        try:
            max_outflow_vel = np.max(disk['Corrected_vel_mag'][outflow_inds])
        except:
            max_outflow_vel = yt.YTQuantity(np.nan, 'km/s')
        sink_dict['max_outflow_speed'].append(max_outflow_vel.in_units('km/s'))
        
        file = open(pickle_file, 'wb')
        pickle.dump((sink_dict), file)
        file.close()
        
        print('updated', pickle_file, 'for file', files.index(fn), 'of', len(files))
    
sys.stdout.flush()
CW.Barrier()

#gather_pickles
if rank == 0:
    pickle_files = sorted(glob.glob('burst_analysis_sink_'+str(sink_id)+'_*.pkl'))

    sink_all = {'time':[], 'mass':[], 'mdot':[], 'max_outflow_speed':[], 'mean_density':[], 'inflow_mass':[], 'outflow_distribution':[] }
    
    for pickle_file in pickle_files:
        file = open(pickle_file, 'rb')
        sink_dict = pickle.load(file)
        file.close()
        
        for key in sink_dict.keys():
            sink_all[key] = sink_all[key] + sink_dict[key]
        
    sorted_inds = np.argsort(sink_all['time'])
    for key in sink_all:
        sink_all[key] = yt.YTArray(sink_all[key])[sorted_inds]
        
    file = open('gathered_burst_analysis_sink_'+str(sink_id)+'.pkl', 'wb')
    pickle.dump((sink_all), file)
    file.close()
    
    print('Gathered sink data for sink_id', sink_id)
    
    file = open('gathered_burst_analysis_sink_'+str(sink_id)+'.pkl', 'rb')
    sink_all = pickle.load(file)
    file.close()
    
    #Maybe plot interesting things?
    smoothing_inds = 5
    mdot_smooth = []
    M_disk_smooth = []
    max_speed_smooth = []
    
    for ind in range(len(sink_all['time'])):
        start_ind = ind - smoothing_inds
        if start_ind < 0:
            if ind == 0:
                mdot_smooth.append(sink_all['mdot'].in_units('msun/yr')[0])
                M_disk_smooth.append(sink_all['inflow_mass'].in_units('msun')[0])
                max_speed_smooth.append(sink_all['max_outflow_speed'].in_units('km/s')[0])
            else:
                mdot_smooth.append(np.nanmean(sink_all['mdot'].in_units('msun/yr')[0:ind]))
                M_disk_smooth.append(np.nanmean(sink_all['inflow_mass'].in_units('msun')[0:ind]))
                max_speed_smooth.append(np.nanmean(sink_all['max_outflow_speed'].in_units('km/s')[0:ind]))
        else:
            mdot_smooth.append(np.nanmean(sink_all['mdot'].in_units('msun/yr')[start_ind:ind]))
            M_disk_smooth.append(np.nanmean(sink_all['inflow_mass'].in_units('msun')[start_ind:ind]))
            max_speed_smooth.append(np.nanmean(sink_all['max_outflow_speed'].in_units('km/s')[start_ind:ind]))
    
    import matplotlib.pyplot as plt
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472
    font_size = 10
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(single_col_width, single_col_width*1.5), sharex=True)
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.0)
    
    axs[0].plot(sink_all['time'], sink_all['mass'].in_units('msun'))
    axs[1].semilogy(sink_all['time'], mdot_smooth)
    axs[2].semilogy(sink_all['time'], M_disk_smooth)
    axs[3].plot(sink_all['time'], max_speed_smooth)
    axs[4].semilogy(sink_all['time'], sink_all['mean_density'].in_units('g/cm**3'))
    #axs[1].semilogy(sink_all['time'], sink_all['mdot'].in_units('msun/yr'))
    #axs[2].semilogy(sink_all['time'], sink_all['inflow_mass'].in_units('msun'))
    #axs[3].plot(sink_all['time'], sink_all['max_outflow_speed'].in_units('km/s'))
    #axs[3].semilogy(sink_all['time'], sink_all['mean_density'].in_units('g/cm**3'))
        
    axs[0].set_ylabel('Mass (Msun)')
    axs[1].set_ylabel('Mdot (Msun/y)')
    #axs[3].set_ylabel('<dens> (g/cm^3)')
    axs[2].set_ylabel('M_d in (Msun)')
    axs[3].set_ylabel('Max speed (km/s)')
    axs[4].set_ylabel('<dens> (g/cm^3)')
    axs[4].set_xlabel('Time (yr)')
    #axs[3].set_xlim([50000, 130000])
    
    #axs[0].set_ylim(bottom=1)
    #axs[1].set_ylim([2.e-6, 2.e-5])
    #axs[3].set_ylim(bottom=1.e4)
    #axs[2].set_ylim([5.e-4, 2.e-2])
    #axs[3].set_ylim([10, 40])
    
    plt.savefig('Sink_id_'+str(sink_id)+'.pdf', bbox_inches='tight', pad_inches=0.02)
        
