#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW
import h5py
import csv
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import my_ramses_module as mym
import my_ramses_fields as myf
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-dt", "--time_step", help="How frequently do you want to measure the outflow?", default=100, type=int)
    parser.add_argument("-disk", "--measure_disks", help="Do you want to measure disk quantities?", type=str, default="False")
    parser.add_argument("-lb", "--lower_bounds", help="lower bounds of box", type=float, default=100.0)
    parser.add_argument("-rad", "--radius", help="What is the radius of the measuring cylindar", type=float, default=1000.0)
    parser.add_argument("-vthres", "--velocity_threshold", help="what velocity threshold do you want to use to outflowing material?", type=float, default=0.0)
    parser.add_argument("-center", "--disk_center", help="do you want two disks centered on the particles or a volume centered on the CoM?", default='CoM', type=str)
    parser.add_argument("-st", "--start_time", help="default is -500yr", default=-500.0, type=float)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-m_all", "--measure_all", help="do you want to measure everything that that is inflowing/outflowing/disk?", default="False")
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def projected_vector(vector, proj_vector):
    """
    Calculates the projection of vecter projected onto vector
    """
    vector_units = vector.units
    proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
    proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
    proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
    proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

#===========================================================================================================================
rank = CW.Get_rank()
size = CW.Get_size()

#define files to read
sys_args = sys.argv
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)
    
args = parse_inputs()
output_file = args.output_file

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = '100'

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')
mym.set_units(units_override)

files = sorted(glob.glob(input_dir+"*/info*.txt"))

#define analysis cylinder:
if args.measure_disks != "False":
    if args.disk_center == 'CoM':
        radius = 100.0
        height = 100.0
        tube_center_1 = None
    else:
        radius = 20.0
        height = 20.0
        tube_center_1 = None
        tube_center_2 = None
else:
    radius = args.radius
    height = 100
    upper_cylinder_bound = args.lower_bounds
    lower_cylinder_bound = -1 * args.lower_bounds
    center_1 = [0.0, 0.0, (upper_cylinder_bound+(height/2.))*yt.units.AU.in_units('cm')]
    center_2 = [0.0, 0.0, (lower_cylinder_bound-(height/2.))*yt.units.AU.in_units('cm')]

#read in sink creation time
ds = yt.load(files[-1], units_override=units_override)
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_centred_sink_id(sink_id)
myf.set_com_pos_use_gas(False)
myf.set_com_vel_use_part(False)
sink_form_time = dd['sink_particle_form_time'][sink_id]
sink_form_companion = dd['sink_particle_form_time'][sink_id+1]#Assumes the companion has the sink id of the primary+1

del dd

if sink_id == 48:
    m_times = mym.generate_frame_times(files, args.time_step, start_time=args.start_time, presink_frames=0, end_time=args.end_time, form_time=sink_form_companion.value)
else:
    m_times = mym.generate_frame_times(files, args.time_step, start_time=args.start_time, presink_frames=0, end_time=args.end_time, form_time=sink_form_time.value)
    
no_frames = len(m_times)
frames = list(range(no_frames))
    
sys.stdout.flush()
CW.Barrier()

verbatim = False
if rank == 0:
    verbatim = True
usable_files = mym.find_files(m_times, files, sink_form_time,sink_id, verbatim=False)
if rank == 0:
    print("Usable_files=", usable_files)
sink_form_time
del files
    
sys.stdout.flush()
CW.Barrier()

#open files to read data out to:
if rank == 0:
    if os.path.isfile(save_dir + output_file) == False:
        f = open(save_dir + output_file, 'a+')
        f.close()
    f = open(save_dir + output_file, 'w')
    if args.measure_disks != "False":
        f.write('Time, Lz_1, Lz_2 \n')
    elif args.measure_all != "False":
        f.write('Time, Particles, Disk, Infalling, Outflowing \n')
    else:
        f.write('Time, Mass, Momentum, Angular Momentum, Max speed, Unbound Mass, CoM dist, Mean speed, F_rad, F_tan \n')
    f.close()

storage = {}
ts = yt.DatasetSeries(usable_files, parallel=True, units_override=units_override)
for sto, ds in ts.piter(storage=storage):
    dd = ds.all_data()
    if args.measure_disks != "False":
        if args.disk_center == 'CoM':
            center_1 = [dd['CoM'][0].in_units('cm'), dd['CoM'][1].in_units('cm'), dd['CoM'][2].in_units('cm')]
            center_2 = None
        else:
            center_1 = [dd['particle_posx'][0].in_units('cm'), dd['particle_posy'][0].in_units('cm'), dd['particle_posz'][0].in_units('cm')]
            center_2 = [dd['particle_posx'][1].in_units('cm'), dd['particle_posy'][1].in_units('cm'), dd['particle_posz'][1].in_units('cm')]
    
    #center_pos = dd['Center_Position'].in_units('AU')
    #center_vel = dd['Center_Velocity'].in_units('cm/s')
    
    part_posx = dd['sink_particle_posx'][sink_id:].in_units('AU')
    part_posy = dd['sink_particle_posy'][sink_id:].in_units('AU')
    part_posz = dd['sink_particle_posz'][sink_id:].in_units('AU')
    part_mass = dd['sink_particle_mass'][sink_id:].in_units('Msun')
    center_pos = yt.YTArray([np.sum(part_posx*part_mass)/np.sum(part_mass), np.sum(part_posy*part_mass)/np.sum(part_mass), np.sum(part_posz*part_mass)/np.sum(part_mass)], 'AU')
    
    part_velx = dd['sink_particle_velx'][sink_id:].in_units('km/s')
    part_vely = dd['sink_particle_vely'][sink_id:].in_units('km/s')
    part_velz = dd['sink_particle_velz'][sink_id:].in_units('km/s')
    center_vel = yt.YTArray([np.sum(part_velx*part_mass)/np.sum(part_mass), np.sum(part_vely*part_mass)/np.sum(part_mass), np.sum(part_velz*part_mass)/np.sum(part_mass)], 'km/s')
    
    rel_pos_x = part_posx - center_pos[0]
    rel_pos_y = part_posy - center_pos[1]
    rel_pos_z = part_posz - center_pos[2]
    rel_vel_x = part_velx - center_vel[0]
    rel_vel_y = part_vely - center_vel[1]
    rel_vel_z = part_velz - center_vel[2]
    
    L_x = part_mass*(rel_pos_y*rel_vel_z - rel_pos_z*rel_vel_y)
    L_y = part_mass*(rel_pos_x*rel_vel_z - rel_pos_z*rel_vel_x)
    L_z = part_mass*(rel_pos_x*rel_vel_y - rel_pos_y*rel_vel_x)
    
    L_orb = yt.YTArray([np.sum(L_x), np.sum(L_y), np.sum(L_z)], 'AU*Msun*km/s')
    L_mag = np.sqrt(np.sum(L_orb[0]**2 + L_orb[1]**2 + L_orb[2]**2))
    L_norm = (L_orb/L_mag).value
    
    #L_orb = yt.YTArray([np.sum(dd['Orbital_Angular_Momentum_x']), np.sum(dd['Orbital_Angular_Momentum_y']), np.sum(dd['Orbital_Angular_Momentum_z'])], 'g*cm**2/s')
    #L_mag = np.sqrt(np.sum(L_orb[0]**2 + L_orb[1]**2 + L_orb[2]**2))
    #L_norm = (L_orb/L_mag).value
    
    center_1 = center_pos + yt.YTArray(L_norm*(height/2. + args.lower_bounds), 'AU')
    center_2 = center_pos + yt.YTArray(-1*L_norm*(height/2 + args.lower_bounds), 'AU')
    
    time_val = ds.current_time.in_units('yr').value - sink_form_time.value

    #define cylinders:
    tube_1 = ds.disk(center_1, L_norm, (radius, 'au'), (height/2., 'au'))
    if center_2 != None:
        tube_2 = ds.disk(center_2, L_norm, (radius, 'au'), (height/2., 'au'))

    #calculate mass in cylinders
    # for region 1
    if args.measure_disks != "False":
        if center_2 != None:
            myf.set_center(1)
            #print "doing center 1, with center_vel =", dd['Center_Velocity']
            Lz_1 = np.sum(tube_1['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
            myf.set_center(2)
            #print "doing center 2, with center_vel =", dd['Center_Velocity']
            Lz_2 = np.sum(tube_2['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
            
            #send data to rank 0 to append to write out.
            write_data = [time_val, Lz_1, Lz_2]
        else:
            myf.set_center(0)
            myf.set_coordinate_system('sph')
            v_kep = tube_1['Relative_Keplerian_Velocity']
            keplerian_inds = np.where((v_kep>0.8) & (v_kep<1.2))[0]
            disk_mass = np.sum(tube_1['cell_mass'][keplerian_inds].in_units('msun'))
            write_data = [time_val, disk_mass]
    elif args.measure_all != "False":
        disk_volume = ds.disk([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], (radius, 'au'), (upper_cylinder_bound, 'au'))
        ang_disk = np.sum(disk_volume['Angular_Momentum'].value)
        if ('all', 'particle_mass') in ds.field_list:
            ang_part = np.sum(disk_volume['Particle_Angular_Momentum'].value)
        else:
            ang_part = 0.0
        myf.set_center(0)
        pos_pos = np.where(tube_1['velz'].in_units('km/s') > args.velocity_threshold)[0]
        other_pos_pos = np.where(tube_1['velz'].in_units('km/s') < args.velocity_threshold)[0]
        ang_in_1 = np.sum(tube_1['Angular_Momentum'][other_pos_pos].value)
        if len(pos_pos) > 0:
            ang_out_1 = np.sum(tube_1['Angular_Momentum'][pos_pos].value)
        else:
            ang_out_1 = 0.0

        neg_pos = np.where(tube_2['velz'].in_units('km/s') < -1*args.velocity_threshold)[0]
        other_neg_pos = np.where(tube_2['velz'].in_units('km/s') > -1*args.velocity_threshold)[0]
        ang_in_2 = np.sum(tube_2['Angular_Momentum'][other_neg_pos].value)
        if len(pos_pos) > 0:
            ang_out_2 = np.sum(tube_2['Angular_Momentum'][neg_pos].value)
        else:
            ang_out_2 = 0.0
                
        ang_in = ang_in_1 + ang_in_2
        ang_out = ang_out_1 + ang_out_2
        write_data = [time_val, ang_part, ang_disk, ang_in, ang_out]
        print("OUTPUT=", time_val, ang_part, ang_disk, ang_in, ang_out)
    else:
        #myf.set_center(0)
        try:
            vel = yt.YTArray([tube_1['x-velocity'].in_units('cm/s'), tube_1['y-velocity'].in_units('cm/s'), tube_1['z-velocity'].in_units('cm/s')], 'cm/s').T
        except:
            vel = yt.YTArray([0, 0, 0], 'cm/s')
        pos_pos = np.where(np.dot(vel, L_norm) > 0)[0]
        #pos_pos = np.where(tube_1['Is_Unbound']==True)[0]
        if len(pos_pos) == 0:
            mass_1 = 0.0
            speed_1 = 0.0
            max_speed_1 = 0.0
            mean_speed_1 = 0.0
            mom_1 = 0.0
            #unbound_1 = 0.0
        else:
            mass_1 = (tube_1['Density'][pos_pos].in_units('msun/au**3')*tube_1['cell_volume'][pos_pos].in_units('au**3')).value
            #mass_1 = tube_1['cell_mass'][pos_pos].in_units('Msun').value
            speed_1 = tube_1['velocity_magnitude'][pos_pos].in_units("km/s").value
            max_speed_1 = np.max(speed_1)
            mean_speed_1 = np.mean(speed_1)
            mom_1 = speed_1*mass_1
            #kin = (tube_1['kinetic_energy'][pos_pos] * tube_1['cell_volume'][pos_pos]).in_units('erg')
            #pot = (tube_1['gpot'][pos_pos]*tube_1['cell_mass'][pos_pos]).in_units('erg')
            #eint = (tube_1['eint'][pos_pos]*tube_1['cell_mass'][pos_pos]).in_units('erg')
            #tot_e = kin + pot# + eint
            #un_ints = np.where(tot_e > 0.0)
            #unbound_1 =tube_1['cell_mass'][pos_pos][un_ints].in_units('Msun').value
        
        #for region 2
        try:
            vel = yt.YTArray([tube_2['x-velocity'].in_units('cm/s'), tube_2['y-velocity'].in_units('cm/s'), tube_2['z-velocity'].in_units('cm/s')], 'cm/s').T
        except:
            vel = yt.YTArray([0, 0, 0], 'cm/s')
        neg_pos = np.where(np.dot(vel, L_norm) < 0)[0]
        #neg_pos = np.where(tube_2['Is_Unbound']==True)[0]
        if len(neg_pos) == 0:
            mass_2 = 0.0
            speed_2 = 0.0
            max_speed_2 = 0.0
            mean_speed_2 = 0.0
            mom_2 = 0.0
            #unbound_2 = 0.0
        else:
            mass_2 = (tube_2['Density'][neg_pos].in_units('msun/au**3')*tube_2['cell_volume'][neg_pos].in_units('au**3')).value
            speed_2 = tube_2['velocity_magnitude'][neg_pos].in_units("km/s").value
            max_speed_2 = np.max(speed_2)
            mean_speed_2 = np.mean(speed_2)
            mom_2 = speed_2*mass_2
            #kin = (tube_2['kinetic_energy'][neg_pos] * tube_2['cell_volume'][neg_pos]).in_units('erg')
            #pot = (tube_2['gpot'][neg_pos]*tube_2['cell_mass'][neg_pos]).in_units('erg')
            #eint = (tube_2['eint'][neg_pos]*tube_2['cell_mass'][neg_pos]).in_units('erg')
            #tot_e = kin + pot# + eint
            #un_ints = np.where(tot_e > 0.0)
            #unbound_2 =tube_2['cell_mass'][neg_pos][un_ints].in_units('Msun').value
        
        #Calculate outflow mass:
        outflow_mass = np.sum(mass_1) + np.sum(mass_2)
        #unbound_mass = np.sum(unbound_1) + np.sum(unbound_2)
        
        #Calculate max outflow speed:
        if max_speed_1 > max_speed_2:
            max_speed = abs(max_speed_1)
        else:
            max_speed = abs(max_speed_2)
        mean_speed = (mean_speed_1 + mean_speed_2)/2.

        #Calculate outflow momentum:
        mom = np.sum(mom_1) + np.sum(mom_2)

        #Calculate outflow angular momentum:
        #L_x = np.sum(tube_1['Angular_Momentum_x'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_x'][neg_pos].in_units("Msun * km**2 / s").value)
        #L_y = np.sum(tube_1['Angular_Momentum_y'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_y'][neg_pos].in_units("Msun * km**2 / s").value)
        #L_z = np.sum(tube_1['Angular_Momentum_z'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_z'][neg_pos].in_units("Msun * km**2 / s").value)
        #L = L_z
        #L = np.sqrt((L_x)**2. + (L_y)**2. + (L_z)**2.)
        
        #Append quantities
        if outflow_mass == 0.0:
            outflow_mass = np.nan
            max_speed = np.nan
            mom = np.nan
            #L = np.nan
            #unbound_mass = np.nan


        print("OUTPUT=", time_val, outflow_mass, mom, max_speed, mean_speed)
            
        #send data to rank 0 to append to write out.
        write_data = [time_val, outflow_mass, mom, max_speed, mean_speed]

    write_string = ''
    for dit in range(len(write_data)):
        if dit != len(write_data)-1:
            write_string = write_string + str(write_data[dit]) + ','
        else:
            write_string = write_string + str(write_data[dit]) + '\n'

    sto.result = write_string
    sto.result_id = str(time_val)
    #del write_data
    #del disk_volume
    #del tube_1
    #del tube_2
    #del pos_pos
    #del other_pos_pos
    #del neg_pos
    #del other_neg_pos
    
#if rank == 0:
pickle_file = save_dir + "AMB.pkl"
file = open(pickle_file, 'wb')
pickle.dump((storage), file)
file.close()

if rank == 0:
    del_keys = []
    for key, value in storage.items():
        if value is None:
            del_keys.append(key)
    for key in del_keys:
        del storage[key]
    f = open(save_dir + output_file, 'a')
    for it in sorted(np.array(list(storage.keys())).astype(np.float)):
        f.write(storage[str(it)])
        print("Printed line:", storage[str(it)])
    f.close()
print("Completed job on rank", rank)
