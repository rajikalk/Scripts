#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields_short as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
cm = plt.cm.get_cmap('seismic')

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("-make_pickles", "--make_pickle_files", type=str, default="True")
    parser.add_argument("-make_plots", "--make_plot_figures", type=str, default="True")
    parser.add_argument("-sphere_radius", "--sphere_radius_cells", type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
        
def projected_vector(vector, proj_vector):
    """
    Calculates the projection of vecter projected onto vector
    """
    vector_units = vector.units
    if len(proj_vector)>3:
        #Calc vector.proj
        v_dot_pv = vector.T[0]*proj_vector.T[0] + vector.T[1]*proj_vector.T[1] + vector.T[2]*proj_vector.T[2]
        pv_dot_pv = proj_vector.T[0]**2 + proj_vector.T[1]**2 + proj_vector.T[2]**2
        proj_v_x = (v_dot_pv/pv_dot_pv)*proj_vector.T[0]
        proj_v_y = (v_dot_pv/pv_dot_pv)*proj_vector.T[1]
        proj_v_z = (v_dot_pv/pv_dot_pv)*proj_vector.T[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T)
    else:
        proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
        proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
        proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

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

if args.make_pickle_files == "True":
    #File files
    usable_files = sorted(glob.glob(input_dir+"*/info*.txt"))

    sys.stdout.flush()
    CW.Barrier()

    #Define units to override:
    units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
    units_override.update({"mass_unit":(2998,"Msun")})
    units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
        
    scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
    scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
    scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
    scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
    mym.set_units(units_override)
    if rank == 0:
        print("set units")

    #find sink particle to center on and formation time
    del units_override['density_unit']
    ds = yt.load(usable_files[-1], units_override=units_override)
    if rank == 0:
        print("Doing initial ds.all_data() load")
    dd = ds.all_data()
    dx_min = np.min(dd['dx'].in_units('au'))
    if args.sphere_radius_cells != None:
        sphere_radius_1 = yt.YTQuantity(args.sphere_radius_cells, 'au')
        sphere_radius_2 = yt.YTQuantity(args.sphere_radius_cells, 'au')
    r_acc = 4*dx_min
    if args.sink_number == None:
        sink_id = np.argmin(dd['sink_particle_speed'])
    else:
        sink_id = args.sink_number
    if rank == 0:
        print("CENTERED SINK ID:", sink_id)
    #myf.set_centred_sink_id(sink_id)
    sink_form_time = dd['sink_particle_form_time'][sink_id]
    del dd
    start_file = mym.find_files([0.0], usable_files, sink_form_time,sink_id)
    usable_files = usable_files[usable_files.index(start_file[0]):]
        
    sys.stdout.flush()
    CW.Barrier()
    
    pickle_file = save_dir + "kepl_mass_" + str(rank) + ".pkl"
    if len(glob.glob(pickle_file[:-4]+"*")) == 0:
        save_dict = {'time':[], 'Kep_mass_primary':[], 'Kep_mass_secondary':[], 'L1':[], 'radius_primary':[], 'mass_profile_primary':[], 'radius_secondary':[], 'mass_profile_secondary':[]}
    else:
        import pdb
        pdb.set_trace()

    file_int = -1
    no_files = len(usable_files)
    for fn in yt.parallel_objects(usable_files, njobs=int(size)):
        frame_name = save_dir + "movie_frame_" + ("%06d" % usable_files.index(fn)) + ".jpg"
        if os.path.isfile(frame_name) == False:
            ds = yt.load(fn, units_override=units_override)
            time_val = ds.current_time.in_units('yr') - sink_form_time
            dd = ds.all_data()
            
            #Get secondary position
            primary_position = yt.YTArray([dd['sink_particle_posx'][sink_id-1], dd['sink_particle_posy'][sink_id-1], dd['sink_particle_posz'][sink_id-1]])
            secondary_position = yt.YTArray([dd['sink_particle_posx'][sink_id], dd['sink_particle_posy'][sink_id], dd['sink_particle_posz'][sink_id]])
            
            primary_velocity = yt.YTArray([dd['sink_particle_velx'][sink_id-1], dd['sink_particle_vely'][sink_id-1], dd['sink_particle_velz'][sink_id-1]])
            secondary_velocity = yt.YTArray([dd['sink_particle_velx'][sink_id], dd['sink_particle_vely'][sink_id], dd['sink_particle_velz'][sink_id]])
            
            primary_mass = dd['sink_particle_mass'][sink_id-1].in_units('g')
            secondary_mass = dd['sink_particle_mass'][sink_id].in_units('g')
            
            if args.sphere_radius_cells == None:
                separation = np.sqrt(np.sum((primary_position - secondary_position)**2))
                #sphere_radius = 0.5 * separation.in_units('au')
                
                sphere_radius_1 = separation.in_units('au')/(np.sqrt(secondary_mass/primary_mass) +1)
                sphere_radius_2 = separation.in_units('au') - sphere_radius_1
            del dd
            
            measuring_sphere_primary = ds.sphere(primary_position.in_units('au'), sphere_radius_1)
            measuring_sphere_secondary = ds.sphere(secondary_position.in_units('au'), sphere_radius_2)
            
            #UPDATE FOR MEASURING KEPLERIAN MASS
            #For primary:
            sph_dx = measuring_sphere_primary['x'].in_units('cm') - primary_position[0].in_units('cm')
            sph_dy = measuring_sphere_primary['y'].in_units('cm') - primary_position[1].in_units('cm')
            sph_dz = measuring_sphere_primary['z'].in_units('cm') - primary_position[2].in_units('cm')
            
            sph_radial_vector = yt.YTArray([sph_dx, sph_dy, sph_dz]).T
            sph_radial_vector_mag = np.sqrt(np.sum(sph_radial_vector**2, axis=1))
            sph_radial_vector_unit = yt.YTArray([sph_dx/sph_radial_vector_mag.value, sph_dy/sph_radial_vector_mag.value, sph_dz/sph_radial_vector_mag.value]).T
            
            sph_dvx = measuring_sphere_primary['velocity_x'].in_units('cm/s') - primary_velocity[0].in_units('cm/s')
            sph_dvy = measuring_sphere_primary['velocity_y'].in_units('cm/s') - primary_velocity[1].in_units('cm/s')
            sph_dvz = measuring_sphere_primary['velocity_z'].in_units('cm/s') - primary_velocity[2].in_units('cm/s')
            
            v_mag = np.sqrt(sph_dvx**2 + sph_dvy**2 + sph_dvz**2)
            sph_velocity_vector = yt.YTArray([sph_dvx, sph_dvy, sph_dvz]).T
            
            shape = np.shape(measuring_sphere_primary['x'])
            radial_vel_vec = yt.YTArray(projected_vector(sph_velocity_vector.in_units('km/s'), sph_radial_vector_unit.in_units('km')).value, 'km/s')
            radial_vel_mag = np.sqrt(np.sum(radial_vel_vec**2, axis=1))
            radial_vel_unit = (radial_vel_vec.T/radial_vel_mag).T
            dot_product = sph_radial_vector_unit.in_units('km').T[0]*radial_vel_unit.T[0] + sph_radial_vector_unit.in_units('km').T[1]*radial_vel_unit.T[1] + sph_radial_vector_unit.in_units('km').T[2]*radial_vel_unit.T[2]
            sign = np.sign(dot_product)
            
            rv_mag = radial_vel_mag*sign
            rv_mag = np.reshape(rv_mag, shape)
            
            tang_vel = np.sqrt(v_mag.in_units('km/s')**2 - rv_mag**2)
            
            #v_kep_pot = np.sqrt(abs(measuring_sphere_primary['Potential'])).in_units('km/s')
            v_kep = np.sqrt((primary_mass*yt.units.gravitational_constant.in_cgs())/sph_radial_vector_mag.in_units('km')).in_units('km/s')
            
            rel_kep = tang_vel/v_kep
            near_kep_inds = np.where((rel_kep>0.8)&(rel_kep<1.2))[0]
            near_kep_mass_primary = np.sum(measuring_sphere_primary['mass'][near_kep_inds].in_units('msun'))
            
            bin_size = 1
            rad_bins_primary = np.arange(0, np.max(sph_radial_vector_mag.in_units('au')).value+bin_size, bin_size)
            bin_centers_primary = (rad_bins_primary[1:] + rad_bins_primary[:-1])/2
            rit = 1
            Mass_profile_primary = []
            while rit < len(rad_bins_primary):
                bin_inds = np.where((sph_radial_vector_mag.in_units('au')>rad_bins_primary[rit-1])&(sph_radial_vector_mag.in_units('au')<rad_bins_primary[rit]))[0]
                try:
                    m_shell = np.mean(measuring_sphere_primary['mass'][bin_inds].in_units('msun'))
                except:
                    m_shell = np.nan
                Mass_profile_primary.append(m_shell)
                rit = rit + 1
            
            save_dict['radius_primary'].append(bin_centers_primary)
            save_dict['mass_profile_primary'].append(Mass_profile_primary)
            
            #Secondary mass
            sph_dx = measuring_sphere_secondary['x'].in_units('cm') - secondary_position[0].in_units('cm')
            sph_dy = measuring_sphere_secondary['y'].in_units('cm') - secondary_position[1].in_units('cm')
            sph_dz = measuring_sphere_secondary['z'].in_units('cm') - secondary_position[2].in_units('cm')
            
            sph_radial_vector = yt.YTArray([sph_dx, sph_dy, sph_dz]).T
            sph_radial_vector_mag = np.sqrt(np.sum(sph_radial_vector**2, axis=1))
            sph_radial_vector_unit = yt.YTArray([sph_dx/sph_radial_vector_mag.value, sph_dy/sph_radial_vector_mag.value, sph_dz/sph_radial_vector_mag.value]).T
            
            sph_dvx = measuring_sphere_secondary['velocity_x'].in_units('cm/s') - secondary_velocity[0].in_units('cm/s')
            sph_dvy = measuring_sphere_secondary['velocity_y'].in_units('cm/s') - secondary_velocity[1].in_units('cm/s')
            sph_dvz = measuring_sphere_secondary['velocity_z'].in_units('cm/s') - secondary_velocity[2].in_units('cm/s')
            
            v_mag = np.sqrt(sph_dvx**2 + sph_dvy**2 + sph_dvz**2)
            sph_velocity_vector = yt.YTArray([sph_dvx, sph_dvy, sph_dvz]).T
            
            shape = np.shape(measuring_sphere_secondary['x'])
            radial_vel_vec = yt.YTArray(projected_vector(sph_velocity_vector.in_units('km/s'), sph_radial_vector_unit.in_units('km')).value, 'km/s')
            radial_vel_mag = np.sqrt(np.sum(radial_vel_vec**2, axis=1))
            radial_vel_unit = (radial_vel_vec.T/radial_vel_mag).T
            dot_product = sph_radial_vector_unit.in_units('km').T[0]*radial_vel_unit.T[0] + sph_radial_vector_unit.in_units('km').T[1]*radial_vel_unit.T[1] + sph_radial_vector_unit.in_units('km').T[2]*radial_vel_unit.T[2]
            sign = np.sign(dot_product)
            
            rv_mag = radial_vel_mag*sign
            rv_mag = np.reshape(rv_mag, shape)
            
            tang_vel = np.sqrt(v_mag.in_units('km/s')**2 - rv_mag**2)
            
            #v_kep_pot = np.sqrt(abs(measuring_sphere_primary['Potential'])).in_units('km/s')
            v_kep = np.sqrt((secondary_mass*yt.units.gravitational_constant.in_cgs())/sph_radial_vector_mag.in_units('km')).in_units('km/s')
            
            rel_kep = tang_vel/v_kep
            near_kep_inds = np.where((rel_kep>0.8)&(rel_kep<1.2))[0]
            near_kep_mass_secondary = np.sum(measuring_sphere_secondary['mass'][near_kep_inds].in_units('msun'))
            
            rad_bins_secondary = np.arange(0, np.max(sph_radial_vector_mag.in_units('au')).value+bin_size, bin_size)
            bin_centers_secondary = (rad_bins_secondary[1:] + rad_bins_secondary[:-1])/2
            rit = 1
            Mass_profile_secondary = []
            while rit < len(rad_bins_secondary):
                bin_inds = np.where((sph_radial_vector_mag.in_units('au')>rad_bins_secondary[rit-1])&(sph_radial_vector_mag.in_units('au')<rad_bins_secondary[rit]))[0]
                try:
                    m_shell = np.mean(measuring_sphere_primary['mass'][bin_inds].in_units('msun'))
                except:
                    m_shell = np.nan
                Mass_profile_secondary.append(m_shell)
                rit = rit + 1
            
            save_dict['radius_secondary'].append(bin_centers_secondary)
            save_dict['mass_profile_secondary'].append(Mass_profile_secondary)
                
            xmax = np.max([bin_centers_primary[-1], bin_centers_secondary[-1]])
            ymin = np.min([np.min(Mass_profile_primary), np.min(Mass_profile_secondary)])
            ymax = np.max([np.max(Mass_profile_primary), np.max(Mass_profile_secondary)])
            plt.clf()
            plt.semilogy(bin_centers_primary, Mass_profile_primary, label="Primary")
            plt.semilogy(bin_centers_secondary, Mass_profile_secondary, label="Secondary")
            plt.xlim([0, xmax])
            plt.ylim([1.e-5, 1.e-9])
            plt.legend(loc='best')
            plt.title("Time:"+str(time_val))
            plt.xlabel("Radius (au)")
            plt.ylabel("Mass (msun)")
            plt.savefig(frame_name)
            
            save_dict['time'].append(time_val)
            save_dict['Kep_mass_primary'].append(near_kep_mass_primary)
            save_dict['Kep_mass_secondary'].append(near_kep_mass_secondary)
            save_dict['L1'].append(sphere_radius_1)
            
            file = open(pickle_file, 'wb')
            pickle.dump((save_dict), file)
            file.close()
            print("wrote file", pickle_file)

sys.stdout.flush()
CW.Barrier()

if args.make_plot_figures == "True":
    
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 10
    
    sink_pickle = "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()

    #radial fraction
    time_arr = []
    Kep_mass_primary = []
    Kep_mass_secondary = []
    L1_arr = []
    bin_centers_primary = []
    Mass_profile_primary = []
    bin_centers_secondary = []
    Mass_profile_secondary = []
    
    pickle_files = sorted(glob.glob(save_dir + "kepl_mass_*.pkl"))
    
    for pick_file in pickle_files:
        file_open = open(pick_file, 'rb')
        save_dict = pickle.load(file_open)
        file_open.close()
        
        time_arr = time_arr + save_dict['time']
        Kep_mass_primary = Kep_mass_primary + save_dict['Kep_mass_primary']
        Kep_mass_secondary = Kep_mass_secondary + save_dict['Kep_mass_secondary']
        L1_arr = L1_arr + save_dict['L1']
        bin_centers_primary = bin_centers_primary + save_dict['radius_primary']
        Mass_profile_primary = Mass_profile_primary + save_dict['mass_profile_primary']
        bin_centers_secondary = bin_centers_secondary + save_dict['radius_secondary']
        Mass_profile_secondary = Mass_profile_secondary + save_dict['mass_profile_secondary']
        
    sort_inds = np.argsort(time_arr)
    time_arr = np.array(time_arr)[sort_inds]
    Kep_mass_primary = np.array(Kep_mass_primary)[sort_inds]
    Kep_mass_secondary = np.array(Kep_mass_secondary)[sort_inds]
    L1_arr = np.array(L1_arr)[sort_inds]
    
    start_ind = np.argmin(abs(particle_data['time'] - time_arr[0]))
    end_ind = np.argmin(abs(particle_data['time'] - time_arr[-1]))

    #Plotting

    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width), sharex=True)
    
    #Plot figure
    ax2 = axs[0].twinx()
    ax2.semilogy(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], color='k')
    if args.sphere_radius_cells != None:
        ax2.axhline(y=args.sphere_radius_cells, ls='--')
    else:
        ax2.semilogy(time_arr, L1_arr, ls='--')
    ax2.set_xlim([particle_data['time'][start_ind], particle_data['time'][end_ind]])
    ax2.set_ylim([np.min(particle_data['separation'][start_ind:end_ind]), np.max(particle_data['separation'][start_ind:end_ind])])
    ax2.set_ylabel('Separation (AU)')
    
    axs[0].semilogy(time_arr,Kep_mass_primary, label="Primary")
    axs[0].semilogy(time_arr,Kep_mass_secondary, label="Secondary")
    axs[0].legend()
    axs[0].set_xlabel('Time (yr)')
    axs[0].set_ylabel('Keplerian mass (M$_\odot$)')

    file_name = save_dir+"Keplerian_mass.pdf"
    plt.savefig(file_name, bbox_inches='tight', dpi=300)
    print("Plotted", file_name)
    
    for time_it in np.arange(len(time_arr)):
        frame_name = save_dir + "movie_frame_" + ("%06d" % time_it) + ".jpg"
        if os.path.isfile(frame_name) == False:
            xmax = np.max([bin_centers_primary[time_it][-1], bin_centers_secondary[time_it][-1]])
            plt.clf()
            plt.semilogy(bin_centers_primary[time_it], Mass_profile_primary[time_it], label="Primary")
            plt.semilogy(bin_centers_secondary[time_it], Mass_profile_secondary[time_it], label="Secondary")
            plt.xlim([0, xmax])
            plt.ylim([1.e-5, 1.e-9])
            plt.legend(loc='best')
            plt.title("Time:"+str(time_arr[time_it]))
            plt.xlabel("Radius (au)")
            plt.ylabel("Mass (msun)")
            plt.savefig(frame_name)
            print("made frame", frame_name)
