#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import my_ramses_module as mym

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("-in_pickle","--input_pickle", default='/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl')
    parser.add_argument("-end_time", "--end_burst_time", type=float)
    parser.add_argument("-make_pickles", "--make_pickle_files", type=str, default="True")
    parser.add_argument("-make_plots", "--make_plot_figures", type=str, default="True")
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

time_bounds = [[5575, 5700], [6570, 6720], [7290, 7365], [7850, 7900]]

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
if args.end_burst_time != None:
    end_time = args.end_burst_time
else:
    event_id = int(input_dir.split('Event_')[-1][0]) - 2
    end_time = time_bounds[event_id][-1]

sys.stdout.flush()
CW.Barrier()

if args.make_pickle_files == "True":
    files = sorted(glob.glob(input_dir+"*/info*.txt"))

    sys.stdout.flush()
    CW.Barrier()

    #Define units to override:
    scale_l = yt.YTQuantity(4, 'pc').in_units('au')
    scale_t = yt.YTQuantity(685706129102738.9, "s").in_units('yr') # 4 pc / 0.18 km/s
    if rank == 0:
        print("set units")

    #find sink particle to center on and formation time
    print("read pickle", args.input_pickle)
    sys.stdout.flush()
    file_open = open(args.input_pickle, 'rb')
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
    file_open.close()
    del particle_data['particle_tag'], particle_data['mass'], particle_data['mdot'], particle_data['separation']
    
    #Get accreted tracer particle IDS
    end_burst_file = mym.find_files([end_time], files, sink_form_time, sink_id, verbatim=False)[0]
    usable_files = files[:files.index(end_burst_file)+1]
    ds = yt.load(end_burst_file)
    dd = ds.all_data()
    min_mass = (-1*(sink_id+1))
    accreted_inds_burst = np.where(dd['particle_mass'] == min_mass)[0]
    accreted_ids_burst = dd['particle_identity'][accreted_inds_burst]
    
    end_sim_file =sorted(glob.glob('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/Level_19/Restart/Level_20_corr_dens_thres/data/output_*/info_*.txt'))[-1]
    ds = yt.load(end_sim_file)
    dd = ds.all_data()
    accreted_inds_all = np.where(dd['particle_mass'] == min_mass)[0]
    accreted_ids_all = dd['particle_identity'][accreted_inds_burst]
    
    accrete_ids_other = yt.YTArray(list(set(accreted_ids_all) - set(accreted_ids_burst)), '')
    not_accreted_ids = yt.YTArray(list(set(dd['particle_identity']) - set(accreted_ids_all)), '')
    
    sys.stdout.flush()
    CW.Barrier()

    file_int = -1
    no_files = len(usable_files)
    for fn in yt.parallel_objects(usable_files, njobs=int(size)):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
        make_pickle = False
        pickle_file = save_dir + "movie_frame_" + ("%06d" % file_int + ".pkl")
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
        if make_pickle:
            ds = yt.load(fn)
            time_val = ds.current_time.value*scale_t - sink_form_time
            dd = ds.all_data()
            
            t_ind = np.argmin(abs(particle_data['time'] - time_val))
            particle_position = particle_data['secondary_position'][t_ind]
            pp_code = particle_position.in_units('pc')/scale_l
            
            accreted_inds_burst = np.in1d(dd['particle_identity'].value, accreted_ids_burst.value).nonzero()[0]
            accrete_inds_other = np.in1d(dd['particle_identity'].value, accrete_ids_other.value).nonzero()[0]
            not_accreted_inds = np.in1d(dd['particle_identity'].value, not_accreted_ids.value).nonzero()[0]
            
            relx = (dd['particle_position_x'][accreted_inds_burst].value - pp_code[0].value)*scale_l
            rely = (dd['particle_position_y'][accreted_inds_burst].value - pp_code[1].value)*scale_l
            relz = (dd['particle_position_z'][accreted_inds_burst].value - pp_code[2].value)*scale_l
            
            burst_positions = [relx, rely, relz]
            
            relx = (dd['particle_position_x'][accrete_inds_other].value - pp_code[0].value)*scale_l
            rely = (dd['particle_position_y'][accrete_inds_other].value - pp_code[1].value)*scale_l
            relz = (dd['particle_position_z'][accrete_inds_other].value - pp_code[2].value)*scale_l
            
            other_positions = [relx, rely, relz]
            
            relx = (dd['particle_position_x'][not_accreted_inds].value - pp_code[0].value)*scale_l
            rely = (dd['particle_position_y'][not_accreted_inds].value - pp_code[1].value)*scale_l
            relz = (dd['particle_position_z'][not_accreted_inds].value - pp_code[2].value)*scale_l
            
            not_accreted_positions = [relx, rely, relz]
            
            write_dict = {'time':time_val, 'burst_positions':burst_positions, 'other_positions':other_positions, 'not_accreted_positions':not_accreted_positions}
            
            file = open(pickle_file, 'wb')
            pickle.dump((write_dict), file)
            file.close()
            print("wrote file", pickle_file, "for file_int", file_int, "of", no_files)

sys.stdout.flush()
CW.Barrier()
'''
if args.make_plot_figures == "True":
    import matplotlib.pyplot as plt
    #plt.rcParams['figure.dpi'] = 300
    from matplotlib.colors import LogNorm
    cm = plt.cm.get_cmap('seismic')
    
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 10
    
    sink_pickle = "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()

    #read pickles and make frames!
    pickle_files = sorted(glob.glob(save_dir + "movie_frame_*.pkl"))
    #Get start and end time

    event_ind = int(input_dir.split('_')[-1][0]) - 2
    time_lim = time_bounds[event_ind]
    start_ind = np.argmin(abs(particle_data['time'] - yt.YTQuantity(time_lim[0], 'yr')))
    end_ind = np.argmin(abs(particle_data['time'] - yt.YTQuantity(time_lim[1], 'yr')))
    
    means_pickle = save_dir+"entrainment.pkl"
    
    if os.path.isfile(means_pickle):
        file = open(means_pickle, 'rb')
        means_dict = pickle.load(file)
        file.close()
        
        try:
            time_arr = means_dict['time']
            acc_arr = means_dict['mdot']
            mean_dens_array = means_dict['density']
            mean_radial_velocity = means_dict['radial_speed']
            mean_radial_momentum = means_dict['radial_momentum']
        except:
            time_arr = means_dict[0]
            acc_arr = means_dict[1]
            mean_dens_array = means_dict[2]
            mean_radial_velocity = means_dict[3]
            mean_radial_momentum = means_dict[4]
    else:
        time_arr = []
        acc_arr = []
        mean_dens_array = []
        mean_radial_velocity = []
        mean_radial_momentum = []
    
    xmin = np.nan
    xmax = np.nan
    ymin = np.nan
    ymax = np.nan
    lin_thresh = np.nan
    rit = -1
    fit = 0
    while fit < len(pickle_files):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            plot_pickle = pickle_files[fit]
            file_name = save_dir + plot_pickle[:-3]+'jpg'
            if os.path.isfile(file_name) == False:
                plt.clf()
                plt.figure(figsize=(12,8))
            
                file = open(plot_pickle, 'rb')
                write_dict = pickle.load(file)
                file.close()
                
                curr_it = np.argmin(abs(particle_data['time'] - write_dict['time']))
                
                time_arr.append(write_dict['time'])
                acc_arr.append(write_dict['mdot'])
                mean_dens_array.append(np.mean(write_dict['density']))
                mean_radial_velocity.append(np.mean(write_dict['radial_speed']))
                mean_radial_momentum.append(np.mean(write_dict['radial_momentum']))
                
                if np.isnan(xmin):
                    xmin = np.min(write_dict['density'].value)
                elif np.min(write_dict['density'].value) < xmin:
                    xmin = np.min(write_dict['density'].value)
                
                if np.isnan(xmax):
                    xmax = np.max(write_dict['density'].value)
                elif np.max(write_dict['density'].value) > xmax:
                    xmax = np.max(write_dict['density'].value)
                    
                if np.isnan(ymin):
                    ymin = np.nanmin(write_dict['radial_momentum'].value)
                elif np.nanmin(write_dict['radial_momentum'].value) < ymin:
                    ymin = np.nanmin(write_dict['radial_momentum'].value)
                    
                if np.isnan(ymax):
                    ymax = np.nanmax(write_dict['radial_momentum'].value)
                elif np.nanmax(write_dict['radial_momentum'].value) > ymax:
                    ymax = np.nanmax(write_dict['radial_momentum'].value)
                    
                if np.isnan(lin_thresh):
                    lin_thresh = np.nanmin(np.abs(write_dict['radial_momentum'].value))
                elif np.nanmin(np.abs(write_dict['radial_momentum'].value)) < lin_thresh:
                    lin_thresh = np.nanmin(np.abs(write_dict['radial_momentum'].value))


                plt.subplot(2,2,1)
                plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind])
                plt.scatter(particle_data['time'][curr_it], particle_data['separation'][curr_it])
                plt.xlim([particle_data['time'][start_ind], particle_data['time'][end_ind]])
                plt.ylim([np.min(particle_data['separation'][start_ind:end_ind]), np.max(particle_data['separation'][start_ind:end_ind])])
                plt.xlabel('Time (yr)')
                plt.ylabel('Separation (au)')
                
                plt.subplot(2,2,3)
                plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'][start_ind:end_ind])
                plt.scatter(particle_data['time'][curr_it], particle_data['mdot'][curr_it][0])
                plt.scatter(particle_data['time'][curr_it], particle_data['mdot'][curr_it][1])
                plt.xlim([particle_data['time'][start_ind], particle_data['time'][end_ind]])
                plt.ylim([np.min(particle_data['mdot'][start_ind:end_ind]), np.max(particle_data['mdot'][start_ind:end_ind])])
                plt.xlabel('Time (yr)')
                plt.ylabel('Separation (au)')
                
                plt.subplot(2,2,(2,4))
                plt.xscale('log')
                plt.yscale('symlog', linthresh=lin_thresh)
                plot = plt.scatter(write_dict['density'], write_dict['radial_momentum'])
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plt.xlabel('density (g/cm$^3$)')
                plt.ylabel('radial momentum (cm$\,$g/s)')

                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                print("Plotted", file_name, "for pickle", fit, "of", len(pickle_files))
                
        fit = fit + 1

means_dict = {'time': time_arr, 'mdot': acc_arr, 'density':mean_dens_array, 'radial_speed':mean_radial_velocity, 'radial_momentum':mean_radial_momentum}
file = open(save_dir+"entrainment.pkl", 'wb')
pickle.dump((means_dict), file)
file.close()

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, 2*single_col_width), sharex=True)
plt.subplots_adjust(hspace=0.0)

axs[0].semilogy(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind])
axs[0].set_xlim(time_bounds[event_ind])
axs[0].set_ylabel('Separation (au)')

axs[1].semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'][start_ind:end_ind])
axs[1].set_ylabel('Accretion (msun/yr)')

axs[2].semilogy(means_dict['time'], means_dict['density'])
axs[2].set_ylabel('Mean_density (g/cm$^3$)')

axs[3].semilogy(means_dict['time'], -1*np.array(means_dict['radial_speed']))
axs[3].set_ylabel('- Mean_radial_velocity (km/s)')

axs[4].semilogy(means_dict['time'], -1*np.array(means_dict['radial_momentum']))
axs[4].set_ylabel('- Mean_radial_momenutm (g*km/s)')
axs[4].set_xlabel('Time (yr)')

plt.savefig('Mean_density.png', bbox_inches='tight')
'''
