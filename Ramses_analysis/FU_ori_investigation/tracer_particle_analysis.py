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
import gc

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

time_bounds = [[3810, 4950], [5575, 5700], [6580, 6730], [7295, 7340], [7850, 7900]]

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
if args.end_burst_time != None:
    end_time = args.end_burst_time
else:
    event_id = int(input_dir.split('Event_')[-1][0]) - 1
    end_time = time_bounds[event_id][-1]

sys.stdout.flush()
CW.Barrier()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
scale_l = yt.YTQuantity(4, 'pc').in_units('au')
scale_t = yt.YTQuantity(685706129102738.9, "s").in_units('yr') # 4 pc / 0.18 km/s
mym.set_units(units_override)
del units_override
gc.collect()

if args.make_pickle_files == "True":
    files = sorted(glob.glob(input_dir+"*/info*.txt"))

    sys.stdout.flush()
    CW.Barrier()

    #Define units to override:
    if rank == 0:
        print("set units")

    #find sink particle to center on and formation time
    if rank == 0:
        print("reading pickle", args.input_pickle)
        sys.stdout.flush()
        file_open = open(args.input_pickle, 'rb')
        particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
        file_open.close()
        del particle_data['particle_tag'], particle_data['mass'], particle_data['mdot'], particle_data['separation'], counter
        gc.collect()
    else:
        sink_id = None
        sink_form_time = None
    sys.stdout.flush()
    CW.Barrier()
    
    if rank != 0:
        sink_id = CW.bcast(sink_id, root=0)
        sink_form_time = CW.bcast(sink_form_time, root=0)
        print("received sink data on rank", rank)
#    try:
#        particle_data = {}
#        particle_data = CW.bcast(particle_data, root=0)
#    except:
#        pass
    

    if rank != 0:
        particle_data = None
        particle_data = CW.bcast(particle_data, root=0)
        print("received particle_data on rank", rank)

    
#    try:
#        particle_data = {}
#        particle_data = CW.Bcast(particle_data, root=0)
#    except:
#        pass
        
#    try:
#        particle_data = None
#        particle_data = CW.Bcast(particle_data, root=0)
#    except:
#        pass
        
    sys.stdout.flush()
    CW.Barrier()
    
    #if rank == 0:
    #Get accreted tracer particle IDS
    print("Getting burst files")
    sys.stdout.flush()
    end_burst_file = mym.find_files([end_time], files, sink_form_time, sink_id, verbatim=True)[0]
    end_file = mym.find_files([end_time+100], files, sink_form_time, sink_id, verbatim=False)[0]
    #end_burst_file = files[-1] #WARNING THIS SHOULD USE THE MYM.FIND FILES LINE
    #end_file = end_burst_file
    print("starting to load end_burs_file")
    ds = yt.load(end_burst_file)
    print("loaded burst file")
    sys.stdout.flush()
    particle_mass = ds.r['particle_mass']
    particle_identity = ds.r['particle_identity']
    print("loaded all data")
    sys.stdout.flush()
    min_mass = (-1*(sink_id+1))
    accreted_inds_burst = np.where(particle_mass == min_mass)[0]
    accreted_ids_burst = particle_identity[accreted_inds_burst]
    print("got burst indexes")
    sys.stdout.flush()
    del particle_mass, particle_identity
    gc.collect()
    
    sys.stdout.flush()
    CW.Barrier()
    #end_sim_file =sorted(glob.glob('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/Level_19/Restart/Level_20_corr_dens_thres/data/output_*/info_*.txt'))[-1]
    print("getting all accreted tracer inds")
    sys.stdout.flush()
    end_sim_file = files[-1]
    usable_files = files[:files.index(end_file)+1]
    ds = yt.load(end_sim_file)
    print("loaded end file")
    sys.stdout.flush()
    particle_mass = ds.r['particle_mass']
    particle_identity = ds.r['particle_identity']

    print("loaded all data")
    sys.stdout.flush()
    accreted_inds_all = np.where(particle_mass == min_mass)[0]
    accreted_ids_all = particle_identity[accreted_inds_all]
    del particle_mass
    gc.collect()
    print("got accreted indexes")
    sys.stdout.flush()
    
    accrete_ids_other = yt.YTArray(list(set(accreted_ids_all.value) - set(accreted_ids_burst.value)), '')
    not_accreted_ids = yt.YTArray(list(set(particle_identity.value) - set(accreted_ids_all.value)), '')
    print('saved other and not accreted tracer particle indices')
    sys.stdout.flush()
    del particle_identity
    gc.collect()
    
    sys.stdout.flush()
    CW.Barrier()

    file_int = -1
    no_files = len(usable_files)
    para_div = 1
    for fn in yt.parallel_objects(usable_files, njobs=int(size/para_div)):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
        make_pickle = False
        pickle_file = save_dir + "movie_frame_" + ("%06d" % file_int + ".pkl")
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
        if make_pickle:
            print('loading file', fn, 'on rank', rank)
            ds = yt.load(fn)
            print('loaded file', fn, 'on rank', rank)
            time_val = ds.current_time.value*scale_t - sink_form_time
            
            t_ind = np.argmin(abs(particle_data['time'] - time_val))
            particle_position = particle_data['secondary_position'][t_ind]
            pp_code = particle_position.in_units('pc')/scale_l
            particle_velocity = particle_data['secondary_velocity'][t_ind]
            pv_code = particle_velocity.in_units('km/s')/scale_v.in_units('km/s')
            
            
            particle_identity = ds.r['gas', 'particle_identity']
            accreted_inds_burst = np.in1d(particle_identity.value, accreted_ids_burst.value).nonzero()[0]
            accrete_inds_other = np.in1d(particle_identity.value, accrete_ids_other.value).nonzero()[0]
            not_accreted_inds = np.in1d(particle_identity.value, not_accreted_ids.value).nonzero()[0]
            del particle_identity
            
            particle_position_x = ds.r['gas', 'particle_position_x']
            particle_position_y = ds.r['gas', 'particle_position_y']
            particle_position_z = ds.r['gas', 'particle_position_z']
            
            
            relx = (particle_position_x[accreted_inds_burst].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[accreted_inds_burst].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[accreted_inds_burst].value - pp_code[2].value)*scale_l

            burst_positions = [relx, rely, relz]
            
            #Get burst velocity
            
            particle_velocity_x = ds.r['gas', 'particle_velocity_x']
            particle_velocity_y = ds.r['gas', 'particle_velocity_y']
            particle_velocity_z = ds.r['gas', 'particle_velocity_z']
            
            rel_vx = (particle_velocity_x[accreted_inds_burst].value - pv_code[0].value)*scale_v.in_units('km/s')
            rel_vy = (particle_velocity_y[accreted_inds_burst].value - pv_code[1].value)*scale_v.in_units('km/s')
            rel_vz = (particle_velocity_z[accreted_inds_burst].value - pv_code[2].value)*scale_v.in_units('km/s')
            del particle_velocity_x, particle_velocity_y, particle_velocity_z
            
            burst_velocity = [rel_vx, rel_vy, rel_vz]
            
            relx = (particle_position_x[accrete_inds_other].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[accrete_inds_other].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[accrete_inds_other].value - pp_code[2].value)*scale_l
            
            other_positions = [relx, rely, relz]
            
            relx = (particle_position_x[not_accreted_inds].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[not_accreted_inds].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[not_accreted_inds].value - pp_code[2].value)*scale_l
            del particle_position_x, particle_position_y, particle_position_z
            
            not_accreted_positions = [relx, rely, relz]
            
            write_dict = {'time':time_val, 'burst_positions':burst_positions, 'other_positions':other_positions, 'not_accreted_positions':not_accreted_positions, 'burst_velocity':burst_velocity}
            
            file = open(pickle_file, 'wb')
            pickle.dump((write_dict), file)
            file.close()
            print("wrote file", pickle_file, "for file_int", file_int, "of", no_files)

sys.stdout.flush()
CW.Barrier()
