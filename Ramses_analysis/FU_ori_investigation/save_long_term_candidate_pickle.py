import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pyramses import rsink
import sys
import os
import yt
import yt.units
from yt.units import g, s, cm, Lsun

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-high_res", "--high_resolution", default=False)
    parser.add_argument("-end_time", "--end_save_time", type=float, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()

top_clean = np.array([177, 292, 48, 51, 262, 195, 17, 10, 75, 159, 272, 275, 176, 118, 54, 45, 85, 103, 71, 101, 258, 150, 93, 221, 151, 154, 102, 168, 175, 56, 309, 239, 109, 73, 72, 83, 141])
top_clean = top_clean[np.argsort(top_clean)]

save_dir = sys.argv[1]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)
    
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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

if args.update_pickle == 'True':
    print("Reading global data")
    if args.high_resolution == 'True':
        file_open = open('/home/100/rlk100/gdata/RAMSES/Global/raw_stars_full_G100_512.pkl', 'rb')
    else:
        file_open = open('/home/100/rlk100/gdata/RAMSES/Global/stars_red_512.pkl', 'rb')
    loaded_sink_data = pickle.load(file_open)
    file_open.close()
    updating = False
    
    sink_ind = args.sink_number
    pickle_name = save_dir+'particle_data_'+str(sink_ind)+'.pkl'
    
    if os.path.isfile(pickle_name):
        try:
            file_open = open(pickle_name, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
        except:
            os.system('cp '+pickle_name.split('.pkl')[0]+'_tmp.pkl '+pickle_name)
            file_open = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating...')
    else:
        updating = True
            
        particle_data = {}
        particle_data.update({'time':np.array([])})
        particle_data.update({'mass':np.array([])})
        particle_data.update({'mdot':np.array([])})
        particle_data.update({'separation':np.array([])})
        particle_data.update({'eccentricity':np.array([])})
        particle_data.update({'closest_sink':np.array([])})
        particle_data.update({'closest_mass':np.array([])})
        particle_data.update({'closest_mdot':np.array([])})
        counter = 0
        sink_form_time = 0
        
    if updating == True:
        loaded_sink_data = loaded_sink_data[counter:]
        for sink_data in loaded_sink_data:
            counter = counter + 1
            if np.remainder(counter, 1000) == 0:
                #try:
                #    os.remove(pickle_name)
                #except:
                #    print("pickle files doesn't exist yet")
                file = open(pickle_name, 'wb')
                pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
                file.close()
                os.system('cp '+pickle_name+' '+pickle_name.split('.pkl')[0]+'_tmp.pkl')
                print('read', counter, 'snapshots of sink particle data, and saved pickle')
            if sink_data[sink_ind][9] > 0:
                if sink_form_time == 0:
                    sink_form_time = sink_data[sink_ind][14]*units['time_unit'].in_units('yr')
                time_val = sink_data[sink_ind][17]*units['time_unit'].in_units('yr') - sink_form_time
                if args.end_save_time == None or time_val < yt.YTQuantity(args.end_save_time, 'yr'):
                    particle_data['time'] = np.append(particle_data['time'], time_val)
                    particle_mass = sink_data[sink_ind][9]*units['mass_unit'].in_units('msun')
                    particle_data['mass'] = np.append(particle_data['mass'], particle_mass)
                    d_mass = sink_data[sink_ind][10]*units['mass_unit'].in_units('msun')
                    d_time = (sink_data[sink_ind][17] - sink_data[sink_ind][18])*units['time_unit'].in_units('yr')
                    acc_val = d_mass/d_time
                    if acc_val == 0:
                        acc_val =1.e-12
                    particle_data['mdot'] = np.append(particle_data['mdot'], acc_val)
                    
                
                    position = np.array([sink_data[sink_ind][0], sink_data[sink_ind][1], sink_data[sink_ind][2]])*units['length_unit'].in_units('au')
                    velocity = np.array([sink_data[sink_ind][3], sink_data[sink_ind][4], sink_data[sink_ind][5]])*units['velocity_unit'].in_units('km/s')
                    closest_ind = np.nan
                    separation = np.inf
                    for other_sink_ind in range(len(sink_data)):
                        if other_sink_ind != sink_ind:
                            if sink_data[other_sink_ind][9]>0:
                                other_pos = np.array([sink_data[other_sink_ind][0], sink_data[other_sink_ind][1], sink_data[other_sink_ind][2]])*units['length_unit'].in_units('au')
                                sep = np.sqrt(np.sum((position - other_pos)**2))
                                if sep < separation:
                                    separation = sep
                                    closest_ind = other_sink_ind
                                
                                
                    particle_data['closest_sink'] = np.append(particle_data['closest_sink'], closest_ind)
                    other_mass = sink_data[closest_ind][9]*units['mass_unit'].in_units('msun')
                    d_mass = sink_data[closest_ind][10]*units['mass_unit'].in_units('msun')
                    d_time = (sink_data[closest_ind][17] - sink_data[closest_ind][18])*units['time_unit'].in_units('yr')
                    acc_val = d_mass/d_time
                    if acc_val == 0:
                        acc_val =1.e-12
                    particle_data['closest_mass'] = np.append(particle_data['closest_mass'], other_mass)
                    particle_data['closest_mdot'] = np.append(particle_data['closest_mdot'], acc_val)
                    particle_data['separation'] = np.append(particle_data['separation'], separation)
                    if np.isnan(closest_ind) == False:
                        other_pos = np.array([sink_data[closest_ind][0], sink_data[closest_ind][1], sink_data[closest_ind][2]])*units['length_unit'].in_units('au')
                        other_vel = np.array([sink_data[closest_ind][3], sink_data[closest_ind][4], sink_data[closest_ind][5]])*units['velocity_unit'].in_units('km/s')
                        CoM_pos = (position*particle_mass + other_pos*other_mass)/(particle_mass + other_mass)
                        CoM_vel = (velocity*particle_mass + other_vel*other_mass)/(particle_mass + other_mass)
                        
                        Cand_vel_rel_to_com = velocity - CoM_vel
                        Other_vel_rel_to_com = other_vel - CoM_vel
                        
                        Cand_pos_rel_to_com = position - CoM_pos
                        Other_pos_rel_to_com = other_pos - CoM_pos
                        
                        Cand_rel_speed_to_com = np.sqrt(np.sum((Cand_vel_rel_to_com)**2))
                        Other_rel_speed_to_com = np.sqrt(np.sum((Other_vel_rel_to_com)**2))
                    
                        reduced_mass = (particle_mass*other_mass)/(particle_mass+other_mass)
                        E_pot = (-1*(yt.units.gravitational_constant_cgs*particle_mass*other_mass)/separation.in_units('cm')).in_units('erg')
                        E_kin = (0.5*particle_mass*Cand_rel_speed_to_com**2 + 0.5*other_mass*Other_rel_speed_to_com**2).in_units('erg')
                        epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                        Cand_r_x_v = yt.YTArray(np.cross(Cand_pos_rel_to_com.in_units('cm'), Cand_vel_rel_to_com.in_units('cm/s')), 'cm**2/s')
                        Other_r_x_v = yt.YTArray(np.cross(Other_pos_rel_to_com.in_units('cm'), Other_vel_rel_to_com.in_units('cm/s')), 'cm**2/s')
                        L = particle_mass.in_units('g')*Cand_r_x_v + other_mass.in_units('g')*Other_r_x_v
                        L_tot = np.sqrt(np.sum(L**2))
                        h_val = L_tot/reduced_mass.in_units('g')
                        e = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.gravitational_constant_cgs*(particle_mass.in_units('g')+other_mass.in_units('g')))**2.))
                    else:
                        e = np.nan
                    particle_data['eccentricity'] = np.append(particle_data['eccentricity'], e)
                
                else:
                    break
        print("Finished saving Sink data")
        #write lastest pickle
        file = open(pickle_name, 'wb')
        pickle.dump((particle_data, counter, sink_ind, sink_form_time), file)
        file.close()
        os.system('cp '+pickle_name+' '+pickle_name.split('.pkl')[0]+'_tmp.pkl')
        print('read', counter, 'snapshots of sink particle data, and saved pickle')
                
'''
Baraffe_mass = yt.YTArray([0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.072, 0.075, 0.080, 0.090, 0.100, 0.110, 0.130, 0.150, 0.170, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400], 'msun')
Baraffe_logL = np.array([-2.469, -2.208, -2.044, -1.783, -1.655, -1.481, -1.399, -1.324, -1.291, -1.261, -1.197, -1.127, -1.154, -1.075, -0.926, -0.795, -0.669, -0.539, -0.199, -0.040, 0.076, 0.171, 0.268, 0.356, 0.436, 0.508, 0.573, 0.634, 0.688, 0.740])
Baraffe_radius = yt.YTArray([0.341, 0.416, 0.472, 0.603, 0.665, 0.796, 0.846, 0.905, 0.942, 0.972, 1.045, 1.113, 1.033, 1.115, 1.270, 1.412, 1.568, 1.731, 2.215, 2.364, 2.458, 2.552, 2.687, 2.821, 2.960, 3.096, 3.227, 3.362, 3.488, 3.621], 'rsun')

#Derive a stellar luminosity
lstar_baraffe_prim = []
rstar_barrafe_prim = []
for mass_val in np.array(particle_data['mass']).T[0]:
    if mass_val < Baraffe_mass[0]:
        lstar_baraffe_prim.append(10**Baraffe_logL[0])
        rstar_barrafe_prim.append(Baraffe_radius[0])
    else:
        closest_inds = sorted(np.argsort(np.abs(Baraffe_mass - mass_val))[:2])
        gradient = (Baraffe_logL[closest_inds][1] - Baraffe_logL[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_logL[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        logL = gradient*mass_val + y_intercept
        lstar_baraffe_prim.append(10**logL)
        
        gradient = (Baraffe_radius[closest_inds][1] - Baraffe_radius[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_radius[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        radius = gradient*mass_val + y_intercept
        rstar_barrafe_prim.append(radius)

lstar_baraffe_prim = yt.YTArray(lstar_baraffe_prim, 'Lsun')
lacc_prim = facc * (yt.units.gravitational_constant_cgs * mass.in_units('g') * mdot.in_units('g/s'))/yt.YTArray(rstar_barrafe_prim).in_units('cm')
ltot_prim = lacc_prim.in_units('lsun') + lstar_baraffe_prim

lstar_baraffe_sec = []
rstar_barrafe_sec = []
for mass_val in particle_data['mass'].T[1]:
    if mass_val < Baraffe_mass[0]:
        lstar_baraffe_sec.append(10**Baraffe_logL[0])
        rstar_barrafe_sec.append(Baraffe_radius[0])
    else:
        closest_inds = sorted(np.argsort(np.abs(Baraffe_mass - mass_val))[:2])
        gradient = (Baraffe_logL[closest_inds][1] - Baraffe_logL[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_logL[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        logL = gradient*mass_val + y_intercept
        lstar_baraffe_sec.append(10**logL)
        
        gradient = (Baraffe_radius[closest_inds][1] - Baraffe_radius[closest_inds][0])/(Baraffe_mass[closest_inds][1] - Baraffe_mass[closest_inds][0])
        y_intercept = Baraffe_radius[closest_inds][1] - gradient*Baraffe_mass[closest_inds][1]
        radius = gradient*mass_val + y_intercept
        rstar_barrafe_sec.append(radius)

lstar_baraffe_sec = yt.YTArray(lstar_baraffe_sec, 'Lsun')
lacc_sec = facc * (yt.units.gravitational_constant_cgs * mass.in_units('g') * mdot.in_units('g/s'))/yt.YTArray(rstar_barrafe_sec).in_units('cm')
ltot_sec = lacc_sec.in_units('lsun') + lstar_baraffe_sec

particle_data.update({'ltot':[ltot_prim, ltot_sec]})

file = open(save_dir+'particle_data_'+str(sink_ind)+'.pkl', 'wb')
pickle.dump((particle_data, counter, sink_ind, sink_form_time), file, protocol=pickle.HIGHEST_PROTOCOL)
file.close()
'''
