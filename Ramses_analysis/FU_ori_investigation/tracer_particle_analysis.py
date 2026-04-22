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
import my_ramses_fields_short as myf
import gc

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("-end_time", "--end_burst_time", type=float)
    parser.add_argument("-make_pickles", "--make_pickle_files", type=str, default="True")
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
    del time_bounds
    gc.collect()

sys.stdout.flush()
CW.Barrier()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
scale_l = yt.YTQuantity(4, 'pc').in_units('au')
scale_t = yt.YTQuantity(685706129102738.9, "s").in_units('yr') # 4 pc / 0.18 km/s

if args.make_pickle_files == "True":
    files = sorted(glob.glob(input_dir+"*/info*.txt"))
    del input_dir
    if os.path.isfile('all_tracer_data.pkl') == False:
        import my_ramses_module as mym
        
        mym.set_units(units_override)

        sys.stdout.flush()
        CW.Barrier()
        
        if args.sink_number != None:
            sink_id = args.sink_number
        else:
            ds = yt.load(files[-1], bbox=bbox)
            sink_id = np.argmin(ds.r['sink_particle_speed'])
        print("SINK ID =", sink_id)
        
        try:
            sink_form_time = ds.r['sink_particle_form_time'][sink_id]
        except:
            bbox= [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]]
            ds = yt.load(files[-1], bbox=bbox)
            sink_form_time = ds.r['sink_particle_form_time'][sink_id]
        
        gc.collect()
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
        print("starting to load end_burst_file")
        ds = yt.load(end_burst_file)
        print("loaded burst file")
        sys.stdout.flush()
        particle_mass = ds.r['particle_mass']
        min_mass = (-1*(sink_id+1))
        accreted_inds_burst = np.where(particle_mass == min_mass)[0]
        del particle_mass
        print("Got all accreted_inds_burst")
        sys.stdout.flush()
        particle_identity = ds.r['particle_identity']
        accreted_ids_burst = particle_identity[accreted_inds_burst]
        print("got burst indexes")
        sys.stdout.flush()
        del particle_identity
        gc.collect()
        
        sys.stdout.flush()
        CW.Barrier()

        print("getting all accreted tracer inds")
        sys.stdout.flush()
        end_sim_file = files[-1]
        usable_files = files[:files.index(end_file)+1]
        ds = yt.load(end_sim_file)
        print("loaded end file")
        sys.stdout.flush()
        particle_mass = ds.r['particle_mass']
        accreted_inds_all = np.where(particle_mass == min_mass)[0]
        del particle_mass
        gc.collect()
        print("got accreted_inds_all")
        sys.stdout.flush()
        particle_identity = ds.r['particle_identity']
        accreted_ids_all = particle_identity[accreted_inds_all]
        print("got accreted indexes")
        sys.stdout.flush()
        
        accrete_ids_other = yt.YTArray(list(set(accreted_ids_all.value) - set(accreted_ids_burst.value)), '')
        not_accreted_ids = yt.YTArray(list(set(particle_identity.value) - set(accreted_ids_all.value)), '')
        print("saved other and not accreted tracer particle indices")
        sys.stdout.flush()
        del particle_identity
        gc.collect()
    
        if rank == 0:
            #Save overall tracer particle data:
            file = open('all_tracer_data.pkl', 'wb')
            pickle.dump((sink_id, sink_form_time, accreted_inds_burst, accreted_ids_burst, accreted_inds_all, accreted_ids_all, accrete_ids_other, not_accreted_ids, end_file), file)
            file.close()
        sys.stdout.flush()
        CW.Barrier()
    else:
        file = open('all_tracer_data.pkl', 'rb')
        sink_id, sink_form_time, accreted_inds_burst, accreted_ids_burst, accreted_inds_all, accreted_ids_all, accrete_ids_other, not_accreted_ids, end_file = pickle.load(file)
        file.close()
        
        usable_files = files[:files.index(end_file)+1]
    
    
    del units_override
    gc.collect()
    sys.stdout.flush()
    CW.Barrier()

    file_int = -1
    no_files = len(usable_files)
    para_div = 1
    for fn in yt.parallel_objects(usable_files, njobs=size):
        print('entering form loop on rank', rank)
        sys.stdout.flush()
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
        make_pickle = False
        pickle_file = save_dir + "movie_frame_" + ("%06d" % file_int + ".pkl")
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
            print("making", pickle_file, "on rank", rank)
        else:
            print(pickle_file, "already exists on rank", rank)
        if make_pickle:
            print('loading file', fn, 'on rank', rank)
            sys.stdout.flush()
            ds = yt.load(fn)
            print('loaded file', fn, 'on rank', rank)
            sys.stdout.flush()
            print("loading tracer particle indices on rank", rank)
            sys.stdout.flush()
            particle_identity = ds.r['particle_identity']
            
            accreted_inds_burst = np.in1d(particle_identity.value, accreted_ids_burst.value).nonzero()[0]
            accrete_inds_other = np.in1d(particle_identity.value, accrete_ids_other.value).nonzero()[0]
            not_accreted_inds = np.in1d(particle_identity.value, not_accreted_ids.value).nonzero()[0]
            del particle_identity
            gc.collect()
        
            print("Getting candidate particle position", rank)
            sys.stdout.flush()
            pp_code = yt.YTArray([ds.r['sink_particle_posx'][sink_id].in_units('code_length'), ds.r['sink_particle_posy'][sink_id].in_units('code_length'), ds.r['sink_particle_posz'][sink_id].in_units('code_length')])
            pv_code = yt.YTArray([ds.r['sink_particle_velx'][sink_id].in_units('code_velocity'), ds.r['sink_particle_vely'][sink_id].in_units('code_velocity'), ds.r['sink_particle_velz'][sink_id].in_units('code_velocity')])
            
                            #Get burst velocity
            print("loading tracer particle velocities on rank", rank)
            sys.stdout.flush()
            curr_accreted_inds = np.where(ds.r['particle_mass'][accreted_inds_burst] == min_mass)[0]
            particle_velocity_x = ds.r['particle_velocity_x'][accreted_inds_burst]
            particle_velocity_x[curr_accreted_inds] = pv_code[0]
            particle_velocity_y = ds.r['particle_velocity_y'][accreted_inds_burst]
            particle_velocity_y[curr_accreted_inds] = pv_code[1]
            particle_velocity_z = ds.r['particle_velocity_z'][accreted_inds_burst]
            particle_velocity_z[curr_accreted_inds] = pv_code[2]
            del curr_accreted_inds
            gc.collect()
            
            rel_vx = (particle_velocity_x.value - pv_code[0].value)*scale_v.in_units('km/s')
            rel_vy = (particle_velocity_y.value - pv_code[1].value)*scale_v.in_units('km/s')
            rel_vz = (particle_velocity_z.value - pv_code[2].value)*scale_v.in_units('km/s')
            del particle_velocity_x, particle_velocity_y, particle_velocity_z
            gc.collect()
            
            burst_velocity = [rel_vx, rel_vy, rel_vz]
            
            print("loading tracer particle positions on rank", rank)
            sys.stdout.flush()
            curr_accreted_inds = np.where(ds.r['particle_mass'] == min_mass)[0]
            particle_position_x = ds.r['particle_position_x']
            particle_position_x[curr_accreted_inds] = pp_code[0]
            particle_position_y = ds.r['particle_position_y']
            particle_position_y[curr_accreted_inds] = pp_code[1]
            particle_position_z = ds.r['particle_position_z']
            particle_position_z[curr_accreted_inds] = pp_code[2]
            
            relx = (particle_position_x[accreted_inds_burst].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[accreted_inds_burst].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[accreted_inds_burst].value - pp_code[2].value)*scale_l

            burst_positions = [relx, rely, relz]
            
            relx = (particle_position_x[accrete_inds_other].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[accrete_inds_other].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[accrete_inds_other].value - pp_code[2].value)*scale_l
            
            other_positions = [relx, rely, relz]
            
            relx = (particle_position_x[not_accreted_inds].value - pp_code[0].value)*scale_l
            rely = (particle_position_y[not_accreted_inds].value - pp_code[1].value)*scale_l
            relz = (particle_position_z[not_accreted_inds].value - pp_code[2].value)*scale_l
            del particle_position_x, particle_position_y, particle_position_z
            gc.collect()
            
            not_accreted_positions = [relx, rely, relz]
            
            time_val = ds.current_time.value*scale_t - sink_form_time
            
            write_dict = {'time':time_val, 'burst_positions':burst_positions, 'other_positions':other_positions, 'not_accreted_positions':not_accreted_positions, 'burst_velocity':burst_velocity}
            
            file = open(pickle_file, 'wb')
            pickle.dump((write_dict), file)
            file.close()
            print("wrote file", pickle_file, "for file_int", file_int, "of", no_files)
            
print("finished saving tracer particle data!")

sys.stdout.flush()
CW.Barrier()
