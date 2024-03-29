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
import numpy.ma as ma
import csv
import re

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-radius", "--proximity_radius", help="within what radius (in AU) do you want to save sink particle data?", type=float, default=10000.0)
    parser.add_argument("-def_sys", "--define_system", help="Is there a particular system that you would like to plot?", type=str, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-end_t", "--end_time", help="dow you want to add a y limit?", type=float, default=None)
    parser.add_argument("-plt_matches", "--plot_matched_times", help="do you want to plot to match times that you found?", default='False', type=str)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def boolean_indexing(v, fillval=np.nan):
    """
    function for squaring off array, so that they are retangular arrays. It filled empty indexes with NaNs by default.
    """
    unit_string = v[0].units
    lens = np.array([len(item) for item in v])
    mask = lens[:,None] > np.arange(lens.max())
    out = np.full(mask.shape,fillval)
    out[mask] = np.concatenate(v)
    out = yt.YTArray(out.T, unit_string)
    return out
    
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]
    
#================================================================================
args = parse_inputs()

path = sys.argv[1]
save_dir = sys.argv[2]
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

if args.update_pickle != 'True':
    print('Not updating pickle')

if args.update_pickle == 'True':
    print("Reading particle data")
    loaded_sink_data = rsink(datadir=path, all=True)
    updating = False

    #Define empty particle array:
    if os.path.isfile('particle_data_raw.pkl'):
        try:
            file_open = open(save_dir+'particle_data_raw.pkl', 'rb')
            particle_data, counter, sink_ind = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            loaded_sink_data = loaded_sink_data[counter:]
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating')
        except:
            os.system('cp '+save_dir+'particle_data_raw_tmp.pkl '+save_dir+'particle_data_raw.pkl ')
            file_open = open(save_dir+'particle_data_raw.pkl', 'rb')
            particle_data, counter, sink_ind = pickle.load(file_open)
            file_open.close()
            counter = int(counter)
            loaded_sink_data = loaded_sink_data[counter:]
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating')
    else:
        updating = True
        if args.sink_number == None:
            last_n = int(sorted(glob.glob(path+"output*"))[-1].split("_")[-1])
            loaded_sink_data_last = rsink(last_n, datadir=path)
            sink_ind = np.argmin(loaded_sink_data_last['u'])
        else:
            sink_ind = args.sink_number
        
        particle_data = {}
        particle_data.update({'particle_tag':[]})
        particle_data.update({'time':[]})
        particle_data.update({'posx':[]})
        particle_data.update({'posy':[]})
        particle_data.update({'posz':[]})
        particle_data.update({'velx':[]})
        particle_data.update({'vely':[]})
        particle_data.update({'velz':[]})
        particle_data.update({'momx':[]})
        particle_data.update({'momy':[]})
        particle_data.update({'momz':[]})
        particle_data.update({'mass':[]})
        particle_data.update({'mdot':[]})
        print("CHECK THAT SEPARATION CALCULATIONS TAKE INTO ACCOUNT CYCLIC BOUNDARIES")
        particle_data.update({'separation':[]})
        particle_data.update({'eccentricity':[]})
        particle_data.update({'E_potential':[]})
        particle_data.update({'E_kinetic':[]})
        particle_data.update({'L_orb':[]})
        particle_data.update({'L_orb_tot':[]})
        particle_data.update({'dist_from_CoM':[]})
        particle_data.update({'period':[]})
        particle_data.update({'semimajor_axis':[]})
        particle_data.update({'semiminor_axis':[]})
        particle_data.update({'last_periastron':[]})
        particle_data.update({'time_since_periastron':[]})
        particle_data.update({'orbital_phase':[]})
        counter = 0
    
    sink_form_time = 0
    last_periastron = np.nan
    found_apastron = True
    prev_separation = 100000
    prev_time = np.nan
    
    if updating == True:
        for sink_data in loaded_sink_data:
            counter = counter + 1
            if np.remainder(counter, 1000) == 0:
                try:
                    os.remove(save_dir+'particle_data_raw.pkl')
                except:
                    print("pickle files doesn't exist yet")
                file = open(save_dir+'particle_data_raw.pkl', 'wb')
                pickle.dump((particle_data, counter, sink_ind), file)
                file.close()
                os.system('cp '+save_dir+'particle_data_raw.pkl '+save_dir+'particle_data_raw_tmp.pkl ')
                print('read', counter, 'snapshots of sink particle data, and saved pickle')
            if len(sink_data['u']) > sink_ind:
                center_pos = np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au')
                dx = sink_data['x']*units['length_unit'].in_units('au') - center_pos[0]
                dy = sink_data['y']*units['length_unit'].in_units('au') - center_pos[1]
                dz = sink_data['z']*units['length_unit'].in_units('au') - center_pos[2]
                dist = np.sqrt(dx**2+dy**2+dz**2)
                particle_tags = np.argwhere(dist < 10000.0).T[0]
                for tag in particle_tags:
                    if tag not in particle_data['particle_tag']:
                        particle_data['particle_tag'].append(tag)
                        print("New particle, tag:", tag)
                mask = []
                for tag in particle_data['particle_tag']:
                    if tag not in particle_tags:
                        mask.append(1)
                    else:
                        mask.append(0)
                masked_tags = ma.masked_array(particle_data['particle_tag'], mask=mask)
                if sink_form_time == 0:
                    sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
                #if len(particle_tags) == 1:
                particle_data['time'].append(time_val)
                particle_data['mass'].append(yt.YTArray(sink_data['m'][masked_tags]*units['mass_unit'].in_units('msun'), 'msun'))
                total_mass = np.sum(particle_data['mass'][-1])
                particle_data['posx'].append(yt.YTArray(sink_data['x'][masked_tags]*units['length_unit'].in_units('au'), 'au'))
                particle_data['posy'].append(yt.YTArray(sink_data['y'][masked_tags]*units['length_unit'].in_units('au'), 'au'))
                particle_data['posz'].append(yt.YTArray(sink_data['z'][masked_tags]*units['length_unit'].in_units('au'), 'au'))
                particle_data['velx'].append(yt.YTArray(sink_data['ux'][masked_tags]*units['velocity_unit'].in_units('km/s'), 'km/s'))
                particle_data['vely'].append(yt.YTArray(sink_data['uy'][masked_tags]*units['velocity_unit'].in_units('km/s'), 'km/s'))
                particle_data['velz'].append(yt.YTArray(sink_data['uz'][masked_tags]*units['velocity_unit'].in_units('km/s'), 'km/s'))
                particle_data['momx'].append(yt.YTArray(sink_data['px'][masked_tags]*units['velocity_unit'].in_units('cm/s')*units['mass_unit'].in_units('g'), 'cm*g/s'))
                particle_data['momy'].append(yt.YTArray(sink_data['py'][masked_tags]*units['velocity_unit'].in_units('cm/s')*units['mass_unit'].in_units('g'), 'cm*g/s'))
                particle_data['momz'].append(yt.YTArray(sink_data['pz'][masked_tags]*units['velocity_unit'].in_units('cm/s')*units['mass_unit'].in_units('g'), 'cm*g/s'))
                
                d_mass = sink_data['dm'][masked_tags]*units['mass_unit'].in_units('msun')
                d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
                particle_data['mdot'].append(yt.YTArray(d_mass/d_time, 'msun/yr'))
                
                position = yt.YTArray([particle_data['posx'][-1],particle_data['posy'][-1], particle_data['posz'][-1]], 'au')
                velocity = yt.YTArray([particle_data['velx'][-1],particle_data['vely'][-1], particle_data['velz'][-1]], 'km/s')
            
                #Calculate center of mass position and velcoity (based on particles only, excluding gas)
                CoM_pos = np.sum((position*particle_data['mass'][-1]).T, axis=0)/np.sum(particle_data['mass'][-1])
                CoM_vel = np.sum((velocity*particle_data['mass'][-1]).T, axis=0)/np.sum(particle_data['mass'][-1])
            
                #Calculate relative positions and velocities and speeds
                pos_rel_to_com = (position.T - CoM_pos).T
                distance_from_com = np.sqrt(np.sum(pos_rel_to_com**2, axis=0))
                particle_data['dist_from_CoM'].append(distance_from_com.in_units('au'))
                
                try:
                    #Calculate relative positions and velocities and speeds
                    vel_rel_to_com = (velocity.T - CoM_vel).T
                        
                    #print("CHECK THAT SEPARATION CALCULATIONS TAKE INTO ACCOUNT CYCLIC BOUNDARIES")
                    separation = np.sum(distance_from_com)
                    if found_apastron == True and separation>prev_separation:
                        last_periastron = prev_time
                        found_apastron = False
                    elif separation<prev_separation and found_apastron == False:
                        found_apastron = True
                    prev_time = time_val
                    prev_separation = separation
                    time_since_p = time_val - last_periastron
                    relative_speed_to_com = np.sqrt(np.sum(vel_rel_to_com**2, axis=0))
                
                    #Calculmate mu and orbital energy
                    reduced_mass = np.product(particle_data['mass'][-1].in_units('g'))/np.sum(particle_data['mass'][-1].in_units('g'))
                    E_pot = (-1*(yt.units.G*np.product(particle_data['mass'][-1].in_units('g')))/separation.in_units('cm')).in_units('erg')
                    E_kin = np.sum((0.5*particle_data['mass'][-1].in_units('g')*relative_speed_to_com.in_units('cm/s')**2).in_units('erg'))
                    
                    primary_ind = np.argmax(particle_data['mass'][-1].in_units('g'))
                    secondary_ind = np.remainder((primary_ind+1),2)
                    
                    #epsilon is the specific orbital energy
                    epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                    
                    #Calculate orbital energy
                    r_x_v = yt.YTArray(np.cross(pos_rel_to_com.T.in_units('cm'), vel_rel_to_com.T.in_units('cm/s')).T, 'cm**2/s')
                    L = particle_data['mass'][-1].in_units('g').T*r_x_v
                    L_tot = np.sqrt(np.sum(np.sum(L, axis=1)**2, axis=0))
                    
                    #h_val is the specific angular momentum
                    h_val = L_tot/reduced_mass.in_units('g')
                    
                    e = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.G*np.sum(particle_data['mass'][-1].in_units('g')))**2.))
                    semimajor_a = ((h_val**2)/(yt.units.G*np.sum(particle_data['mass'][-1].in_units('g'))*(1-e**2))).in_units('AU')
                    period = (2*np.pi*np.sqrt((semimajor_a.in_units('AU')**3)/(yt.units.G*np.sum(particle_data['mass'][-1].in_units('g'))))).in_units('yr')
                    phase = time_since_p/period
                    semiminor_a = semimajor_a * np.sqrt((1-e**2))
                    
                    particle_data['eccentricity'].append(e)
                    particle_data['separation'].append(separation.in_units('au'))
                    particle_data['L_orb'].append(L)
                    particle_data['L_orb_tot'].append(L_tot)
                    particle_data['E_potential'].append(E_pot)
                    particle_data['E_kinetic'].append(E_kin)
                    particle_data['period'].append(period)
                    particle_data['semimajor_axis'].append(semimajor_a)
                    particle_data['semiminor_axis'].append(semiminor_a)
                    particle_data['last_periastron'].append(last_periastron)
                    particle_data['time_since_periastron'].append(time_since_p)
                    particle_data['orbital_phase'].append(phase)
                    
                except:
                    particle_data['eccentricity'].append(yt.YTArray(np.nan, ''))
                    particle_data['separation'].append(yt.YTArray(np.nan, 'au'))
                    particle_data['L_orb'].append(yt.YTArray([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]], 'cm**2*g/s'))
                    particle_data['L_orb_tot'].append(yt.YTArray([np.nan, np.nan, np.nan], 'cm**2*g/s'))
                    particle_data['E_potential'].append(yt.YTArray(np.nan, 'erg'))
                    particle_data['E_kinetic'].append(yt.YTArray(np.nan, 'erg'))
                    particle_data['period'].append(yt.YTArray(np.nan, 'yr'))
                    particle_data['semimajor_axis'].append(yt.YTArray(np.nan, 'au'))
                    particle_data['semiminor_axis'].append(yt.YTArray(np.nan, 'au'))
                    particle_data['last_periastron'].append(yt.YTArray(np.nan, 'yr'))
                    particle_data['time_since_periastron'].append(yt.YTArray(np.nan, 'yr'))
                    particle_data['orbital_phase'].append(yt.YTArray(np.nan, ''))

    print("Finished reading particle data")


try:
    print('Reading neat pickle file')
    file_open = open(save_dir+'particle_data_neat.pkl', 'rb')
    particle_data, counter, sink_ind = pickle.load(file_open)
    file_open.close()
except:
    print('Reading raw pickle file')
    file_open = open(save_dir+'particle_data_raw.pkl', 'rb')
    particle_data, counter, sink_ind = pickle.load(file_open)
    file_open.close()
    
    print("making arrays neat")

    for key in list(particle_data.keys()):
        try:
            particle_data[key] = boolean_indexing(particle_data[key])
        except:
            try:
                unit_string = particle_data[key][0].units
                print('unit_string =', unit_string)
            except:
                unit_string = ""
                print('unit_string =', unit_string)
                
            particle_data[key] = yt.YTArray(np.array(particle_data[key]), unit_string)
            print("1D array can't be squared off")
        
    file = open(save_dir+'particle_data_neat.pkl', 'wb')
    pickle.dump((particle_data, counter, sink_ind), file)
    file.close()
#Save neat and particle_data_neat
'''
plt.clf()
fig = plt.figure()
fig.set_size_inches(6, 10.)
gs = gridspec.GridSpec(4, 1)
gs.update(hspace=0.0)
ax1 = fig.add_subplot(gs[0,0]) #Separation
ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
ax4 = fig.add_subplot(gs[3,0], sharex=ax1)


ax1.plot(particle_data['time'], particle_data['orbital_phase'])
ax1.set_ylabel('Orbital Phase')
plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
ax1.set_xlim([particle_data['time'][0].value, particle_data['time'][-1].value])

ax2.semilogy(particle_data['time'], particle_data['eccentricity'])
ax2.set_ylabel('Eccentricity')
plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)

ax3.semilogy(particle_data['time'], particle_data['period'])
ax3.set_ylabel('Period (yr)')

ax4.semilogy(particle_data['time'], particle_data['separation'])
ax4.set_ylabel('Separation (AU)')
ax4.set_xlabel('Time (yr)')

plt.savefig("paper_data.pdf", bbox_inches='tight', pad_inches=0.02)
print("Created image paper_data.pdf")
'''
image_name = save_dir + "total_system_evolution"
#line_style = [':', '-.', '--', '-']
#line_colour = ['b', 'orange', 'g', 'm']
labels = []
for tag in particle_data['particle_tag']:
    labels.append(r'Sink ' + str(int(tag.value)))

#Plot binary evolution
plt.clf()
fig = plt.figure()
fig.set_size_inches(6, 10.)
gs = gridspec.GridSpec(4, 1)
gs.update(hspace=0.0)
ax1 = fig.add_subplot(gs[0,0]) #Separation
ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #Eccentricity
#ax4 = fig.add_subplot(gs[3,0], sharex=ax1)

for nit in range(len(labels)):
    try:
        ax1.semilogy(particle_data['time'].value, particle_data['dist_from_CoM'].T[nit].value, label=labels[nit], lw=0.5)#, color=line_colour[nit])
    except:
        ax1.semilogy(particle_data['time'].value, particle_data['dist_from_CoM'][nit].value, label=labels[nit], lw=0.5)#, color=line_colour[nit]
ax1.set_ylabel('Distance from CoM (AU)')
ax1.legend(loc='upper left')
ax1.set_ylim(bottom=1.e0)
plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
ax1.set_xlim([particle_data['time'][0].value, particle_data['time'][-1].value])

#Plot accretion
for pit in range(len(labels)):
    try:
        ax2.semilogy(particle_data['time'].value, particle_data['mdot'].value.T[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
    except:
        ax2.semilogy(particle_data['time'].value, particle_data['mdot'].value[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
try:
    total_accretion = np.sum(np.nan_to_num(particle_data['mdot'].value.T), axis=0)
    ax2.semilogy(particle_data['time'].value, total_accretion, color='k', label="Total", lw=0.5)
except:
    total_accretion = np.sum(np.nan_to_num(particle_data['mdot'].value), axis=0)
    ax2.semilogy(particle_data['time'].value, total_accretion, color='k', label="Total", lw=0.5)
ax2.set_ylabel(r'Accretion Rate (M$_\odot$/yr)')
#ax2.legend(loc='best')
plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)

for pit in range(len(labels)):
    try:
        ax3.semilogy(particle_data['time'].value, particle_data['mass'].value.T[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
    except:
        ax3.semilogy(particle_data['time'].value, particle_data['mass'].value[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
try:
    total_mass = np.sum(np.nan_to_num(particle_data['mass'].value.T), axis=0)
    ax3.semilogy(particle_data['time'].value, total_mass, color='k', label="Total", lw=0.5)
except:
    total_mass = np.sum(np.nan_to_num(particle_data['mass'].value), axis=0)
    ax3.semilogy(particle_data['time'].value, total_mass, color='k', label="Total", lw=0.5)
ax3.set_ylabel(r'Mass (M$_\odot$)')
ax3.set_ylim(bottom=1.e-2)
#ax2.legend(loc='best')

#Save image
plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
print("Created image", image_name)


#Anaylsys particular systems
if args.define_system != None:
    systems_hierarchy = args.define_system
else:
    print("Not making plots for the multiple star system")

reduced_systems_data = {}
reduced_systems_data.update({'tag':[]})
reduced_systems_data.update({'base_tags':[]})
reduced_systems_data.update({'posx':[]})
reduced_systems_data.update({'posy':[]})
reduced_systems_data.update({'posz':[]})
reduced_systems_data.update({'velx':[]})
reduced_systems_data.update({'vely':[]})
reduced_systems_data.update({'velz':[]})
reduced_systems_data.update({'mass':[]})
reduced_systems_data.update({'mdot':[]})
reduced_systems_data.update({'mdot_individual':[]})
print("CHECK THAT SEPARATION CALCULATIONS TAKE INTO ACCOUNT CYCLIC BOUNDARIES")
reduced_systems_data.update({'separation':[]})
reduced_systems_data.update({'eccentricity':[]})
reduced_systems_data.update({'E_potential':[]})
reduced_systems_data.update({'E_kinetic':[]})
reduced_systems_data.update({'L_orb_tot':[]})
reduced_systems_data.update({'dist_from_CoM':[]})
reduced_systems_data.update({'period':[]})
reduced_systems_data.update({'semimajor_axis':[]})
reduced_systems_data.update({'semiminor_axis':[]})
reduced_systems_data.update({'last_periastron':[]})
reduced_systems_data.update({'time_since_periastron':[]})
reduced_systems_data.update({'orbital_phase':[]})
reduced_systems_data.update({'time':particle_data['time']})

system_tag_int = 65
if systems_hierarchy != None:
    while len(find(systems_hierarchy, '[')) > 0:
        indexes = []
        for ind, char in enumerate(systems_hierarchy):
            if char == '[':
                indexes.append(ind)
            if char == ']':
                start_ind = indexes.pop()
                end_ind = ind
                
                binary_tags = eval('['+systems_hierarchy[start_ind+1:end_ind]+']')
                system_tag = chr(system_tag_int)
                system_tag_int = system_tag_int + 1
                    
                #define array like position, velocity, mass
                
                binary_data = {}
                binary_data.update({'posx':[]})
                binary_data.update({'posy':[]})
                binary_data.update({'posz':[]})
                binary_data.update({'velx':[]})
                binary_data.update({'vely':[]})
                binary_data.update({'velz':[]})
                binary_data.update({'mass':[]})
                binary_data.update({'mdot':[]})
                
                for tag in binary_tags:
                    if tag in reduced_systems_data['tag']:
                        system_ind = reduced_systems_data['tag'].index(tag)
                        binary_data['mass'].append(reduced_systems_data['mass'][system_ind])
                        binary_data['mdot'].append(reduced_systems_data['mdot'][system_ind])
                        binary_data['posx'].append(reduced_systems_data['posx'][system_ind])
                        binary_data['posy'].append(reduced_systems_data['posy'][system_ind])
                        binary_data['posz'].append(reduced_systems_data['posz'][system_ind])
                        binary_data['velx'].append(reduced_systems_data['velx'][system_ind])
                        binary_data['vely'].append(reduced_systems_data['vely'][system_ind])
                        binary_data['velz'].append(reduced_systems_data['velz'][system_ind])
                    else:
                        system_ind = np.argwhere(particle_data['particle_tag'] == tag)[0][0]
                        binary_data['mass'].append(particle_data['mass'][system_ind])
                        binary_data['mdot'].append(particle_data['mdot'][system_ind])
                        binary_data['posx'].append(particle_data['posx'][system_ind])
                        binary_data['posy'].append(particle_data['posy'][system_ind])
                        binary_data['posz'].append(particle_data['posz'][system_ind])
                        binary_data['velx'].append(particle_data['velx'][system_ind])
                        binary_data['vely'].append(particle_data['vely'][system_ind])
                        binary_data['velz'].append(particle_data['velz'][system_ind])
                        
                for key in list(binary_data.keys()):
                    try:
                        unit_string = binary_data[key][0].units
                        print('unit_string =', unit_string)
                    except:
                        unit_string = ""
                        print('unit_string =', unit_string)
                        
                    binary_data[key] = yt.YTArray(np.array(binary_data[key]), unit_string)
                    
                
                binary_masses = binary_data['mass']
                binary_accretion = binary_data['mdot']
                binary_positions = yt.YTArray([binary_data['posx'].in_units('au'), binary_data['posy'].in_units('au'), binary_data['posz'].in_units('au')], 'au')
                binary_velocities = yt.YTArray([binary_data['velx'].in_units('km/s'), binary_data['vely'].in_units('km/s'), binary_data['velz'].in_units('km/s')], 'km/s')
                
                CoM_pos = np.sum(binary_masses*binary_positions, axis=1)/np.sum(binary_masses, axis=0)
                CoM_vel = np.sum(binary_masses*binary_velocities, axis=1)/np.sum(binary_masses, axis=0)
                
                pos_rel_to_com_1 = binary_positions[:,0] - CoM_pos
                pos_rel_to_com_2 = binary_positions[:,1] - CoM_pos
                distance_from_com_1 = np.sqrt(np.sum(pos_rel_to_com_1**2, axis=0))
                distance_from_com_2 = np.sqrt(np.sum(pos_rel_to_com_2**2, axis=0))
                
                vel_rel_to_com_1 = binary_velocities[:,0] - CoM_vel
                vel_rel_to_com_2 = binary_velocities[:,1] - CoM_vel

                print("CHECK THAT SEPARATION CALCULATIONS TAKE INTO ACCOUNT CYCLIC BOUNDARIES")
                separation = distance_from_com_1 + distance_from_com_2
                relative_speed_to_com_1 = np.sqrt(np.sum(vel_rel_to_com_1**2, axis=0))
                relative_speed_to_com_2 = np.sqrt(np.sum(vel_rel_to_com_2**2, axis=0))
            
                #Calculmate mu and orbital energy
                reduced_mass = np.product(binary_masses, axis=0).in_units('g**2')/np.sum(binary_masses, axis=0).in_units('g')
                E_pot = (-1*(yt.units.G*np.product(binary_masses, axis=0).in_units('g**2'))/separation.in_units('cm')).in_units('erg')
                E_kin = (0.5*binary_masses[0].in_units('g')*relative_speed_to_com_1.in_units('cm/s')**2).in_units('erg') + (0.5*binary_masses[1].in_units('g')*relative_speed_to_com_2.in_units('cm/s')**2).in_units('erg')
                
                #epsilon is the specific orbital energy
                epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                
                #Calculate orbital energy
                r_x_v_1 = yt.YTArray(np.cross(pos_rel_to_com_1.T.in_units('cm'), vel_rel_to_com_1.T.in_units('cm/s')), 'cm**2/s')
                r_x_v_2 = yt.YTArray(np.cross(pos_rel_to_com_2.T.in_units('cm'), vel_rel_to_com_2.T.in_units('cm/s')), 'cm**2/s')
                L_1 = r_x_v_1.T*binary_masses[0].in_units('g')
                L_2 = r_x_v_2.T*binary_masses[1].in_units('g')
                L_tot = np.sqrt(np.sum(L_1**2, axis=0)) + np.sqrt(np.sum(L_2**2, axis=0))
                
                #h_val is the specific angular momentum
                h_val = L_tot/reduced_mass.in_units('g')
                
                e = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.G*np.sum(binary_masses, axis=0).in_units('g'))**2.))
                
                semimajor_a = ((h_val**2)/(yt.units.G*np.sum(binary_data['mass'].in_units('g'),axis=0)*(1-e**2))).in_units('AU')
                period = (2*np.pi*np.sqrt((semimajor_a.in_units('AU')**3)/(yt.units.G*np.sum(binary_data['mass'].in_units('g'),axis=0)))).in_units('yr')
                
                towards_apastron = np.argwhere(attractor>0)
                attractor = np.sign(separation[1:] - separation[:-1])
                away_from_apastron = np.argwhere(attractor<0)
                phase = np.arccos(((((semimajor_a*(1-e**2))/separation) - 1)/e))/np.pi/2
                phase[away_from_apastron] = 1- phase[away_from_apastron]
                
                if 'last_periastron' not in particle_data.keys():
                    particle_data.update({'last_periastron':[]})
                    last_periastron = np.nan
                    found_apastron = True
                    prev_separation = 100000
                    prev_time = particle_data['time'][0]
                    
                    for sep_it in range(len(separation)):
                        if found_apastron == True and separation[sep_it]>prev_separation:
                            last_periastron = prev_time
                            print('Found  Periastron at', prev_time)
                            found_apastron = False
                        elif separation[sep_it]<prev_separation and found_apastron == False:
                            found_apastron = True
                        prev_time = particle_data['time'][sep_it]
                        prev_separation = separation[sep_it]

                        particle_data['last_periastron'].append(last_periastron)
                    
                    '''
                    for sep_it in range(len(particle_data['separation'])):
                        if found_apastron == True and particle_data['separation'][sep_it]>prev_separation:
                            last_periastron = prev_time
                            print('Found  Periastron at', prev_time)
                            found_apastron = False
                        elif particle_data['separation'][sep_it]<prev_separation and found_apastron == False:
                            found_apastron = True
                        prev_time = particle_data['time'][sep_it]
                        prev_separation = particle_data['separation'][sep_it]

                        particle_data['last_periastron'].append(last_periastron)
                    '''
                
                try:
                    time_since_p = particle_data['time'].value - particle_data['last_periastron'].value
                except:
                    time_since_p = particle_data['time'].value - particle_data['last_periastron']
                phase = time_since_p/period.value
                semiminor_a = semimajor_a * np.sqrt((1-e**2))
                
                
                reduced_systems_data['tag'].append(system_tag)
                reduced_systems_data['base_tags'].append(binary_tags)
                reduced_systems_data['posx'].append(CoM_pos[0])
                reduced_systems_data['posy'].append(CoM_pos[1])
                reduced_systems_data['posz'].append(CoM_pos[2])
                reduced_systems_data['velx'].append(CoM_vel[0])
                reduced_systems_data['vely'].append(CoM_vel[1])
                reduced_systems_data['velz'].append(CoM_vel[2])
                reduced_systems_data['mass'].append(np.sum(binary_masses, axis=0))
                reduced_systems_data['mdot'].append(np.sum(binary_accretion, axis=0))
                reduced_systems_data['mdot_individual'].append(binary_accretion)
                print("CHECK THAT SEPARATION CALCULATIONS TAKE INTO ACCOUNT CYCLIC BOUNDARIES")
                reduced_systems_data['separation'].append(separation)
                reduced_systems_data['eccentricity'].append(e)
                reduced_systems_data['E_potential'].append(E_pot)
                reduced_systems_data['E_kinetic'].append(E_kin)
                reduced_systems_data['L_orb_tot'].append(L_tot)
                reduced_systems_data['dist_from_CoM'].append([distance_from_com_1, distance_from_com_2])
                reduced_systems_data['period'].append(period)
                reduced_systems_data['semimajor_axis'].append(semimajor_a)
                reduced_systems_data['semiminor_axis'].append(semiminor_a)
                reduced_systems_data['time_since_periastron'].append(time_since_p)
                reduced_systems_data['orbital_phase'].append(phase)

                systems_hierarchy = systems_hierarchy[:start_ind] + "\'" + system_tag + "\'" + systems_hierarchy[end_ind+1:]
                print('systems_hierarchy =', systems_hierarchy)
                
                file = open(save_dir+'reduced_system_data.pkl', 'wb')
                pickle.dump((reduced_systems_data), file)
                file.close()
                print("wrote file", save_dir+'reduced_system_data.pkl')
                
                break
                
                '''
                if args.plot_matched_times != 'False':
                    plt.clf()
                    fig = plt.figure()
                    matches_dict = {}
                    with open('./quadrupole_rv_times.csv','rU') as match_times:
                        for line in match_times:
                            row = re.split('[[]|[]]', line)
                            if row[0][0] != '#':
                                no_thres = eval('['+row[1]+']')
                                obs_thres = eval('['+row[3]+']')
                                matches_dict.update({str(row[0][0]):[no_thres, obs_thres]})
                    #proj_colours = plt.cm.nipy_spectral(np.linspace(0,1,8))
                    normalize = mcolors.Normalize(vmin=0, vmax=8)
                    colormap =cm.nipy_spectral
                    #proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']
                    plt.semilogy(particle_data['time'].value, separation.value, lw=0.5)
                    plt.xlim(left=particle_data['time'].value[np.where(np.isnan(separation.value) == False)[0][0]])
                    if args.end_time != None:
                        plt.xlim(right=float(args.end_time))
                    for proj in range(8):
                        alpha = [0.25, 1.0]
                        for thres_it in range(2):
                            for match_times in matches_dict[str(proj)][thres_it]:
                                try:
                                    plt.axvspan(match_times[0], match_times[1], alpha=thres_it, color=colormap(normalize(proj)))# color=proj_colours[proj])
                                    print("plotted time range")
                                except:
                                    plt.axvline(match_times, alpha=thres_it, color=colormap(normalize(proj)), lw=2)
                    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
                    scalarmappaple.set_array(np.arange(8))
                    plt.colorbar(scalarmappaple, ticks=np.arange(8), pad=0.0)
                    image_name = save_dir + "match_time_" + str(binary_tags[0]) + "_" + str(binary_tags[1])
                    plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
                    print("Created image", image_name)
                    
                '''
                
                #Plot figures
                image_name = save_dir + "binary_evolution_plot_" + str(binary_tags[0]) + "_" + str(binary_tags[1])
                #line_style = [':', '-.', '--', '-']
                #line_colour = ['b', 'orange', 'g', 'm']
                labels = []
                for tag in binary_tags:
                    if type(tag) == str:
                        system_list = [tag]
                        system_string = str(system_list)
                        while "'" in system_string:
                            s_ind = system_string.index("'")
                            sys_tag = system_string[s_ind+1:s_ind+2]
                            sys_int = reduced_systems_data['tag'].index(sys_tag)
                            base_tags = reduced_systems_data['base_tags'][sys_int]
                            system_string = system_string[:s_ind] + str(base_tags) + system_string[s_ind+3:]
                        label_string = tag + ' ('+system_string[1:-1]+')'
                        labels.append(label_string)
                    else:
                        labels.append(r'Sink ' + str(tag))

                non_nan_inds = np.argwhere(np.isnan(separation) == False).T[0]
    
                #Plot binary evolution
                plt.clf()
                fig = plt.figure()
                fig.set_size_inches(6, 10.)
                #gs = gridspec.GridSpec(4, 1)
                gs = gridspec.GridSpec(3, 1)
                gs.update(hspace=0.0)
                ax1 = fig.add_subplot(gs[0,0]) #Separation
                ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
                ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #Eccentricity
                #ax4 = fig.add_subplot(gs[3,0], sharex=ax1)

                ax1.semilogy(particle_data['time'].value, separation.value, lw=0.5)#, color=line_colour[nit])
                ax1.axhline(y=25, color='k')
                #ax1.set_ylim(bottom=1.e-1)
                plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
                if args.end_time != None:
                    ax1.set_xlim(right=float(args.end_time))
                else:
                    ax1.set_xlim([particle_data['time'][non_nan_inds][0].value, particle_data['time'][non_nan_inds][-1].value])
                ax1.yaxis.set_ticks_position('both')
                ax1.tick_params(axis='y', which='major', direction="in")
                ax1.tick_params(axis='y', which='minor', direction="in")
                ax1.tick_params(axis='x', which='major', direction="in")
                ax1.set_ylabel('Separation (AU)')

                #Plot accretion
                for pit in range(len(labels)):
                    ax2.semilogy(particle_data['time'].value, binary_accretion[pit].value, label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
                total_accretion = np.sum(binary_accretion, axis=0)
                ax2.legend(loc='best')
                ax2.semilogy(particle_data['time'].value, total_accretion.value, color='k', label="Total", lw=0.5)
                ax2.set_ylabel(r'Accretion Rate (M$_\odot$/yr)')
                ax2.yaxis.set_ticks_position('both')
                ax2.tick_params(axis='y', which='major', direction="in")
                ax2.tick_params(axis='y', which='minor', direction="in")
                ax2.tick_params(axis='x', which='major', direction="in")
                #ax2.legend(loc='best')
                plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
                '''
                for pit in range(len(labels)):
                    ax3.semilogy(particle_data['time'].value, binary_masses[pit].value, label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
                total_mass = np.sum(binary_masses, axis=0)
                ax3.semilogy(particle_data['time'].value, total_mass.value, color='k', label="Total", lw=0.5)
                ax3.set_ylabel(r'Mass (M$_\odot$)')
                ax3.set_ylim(bottom=1.e-2)
                ax3.yaxis.set_ticks_position('both')
                ax3.tick_params(axis='y', which='major', direction="in")
                ax3.tick_params(axis='y', which='minor', direction="in")
                ax3.tick_params(axis='x', which='major', direction="in")
                plt.setp([ax3.get_xticklabels() for ax3 in fig.axes[:-1]], visible=False)
                
                ax4.semilogy(particle_data['time'].value, e.value, lw=0.5)
                ax4.set_ylabel(r'Eccentricity')
                #ax4.set_ylim([0.0, 2.0])
                ax4.set_xlabel(r'Time (yr)')
                ax4.yaxis.set_ticks_position('both')
                ax4.tick_params(axis='y', which='major', direction="in")
                ax4.tick_params(axis='y', which='minor', direction="in")
                ax4.tick_params(axis='x', which='major', direction="in")
                '''
                #Save image
                
                ax3.semilogy(particle_data['time'].value, e.value, lw=0.5)
                ax3.set_ylabel(r'Eccentricity')
                #ax4.set_ylim([0.0, 2.0])
                ax3.set_xlabel(r'Time (yr)')
                ax3.yaxis.set_ticks_position('both')
                ax3.tick_params(axis='y', which='major', direction="in")
                ax3.tick_params(axis='y', which='minor', direction="in")
                ax3.tick_params(axis='x', which='major', direction="in")
                
                plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
                print("Created image", image_name)
                
                import pdb
                pdb.set_trace()
                
                image_name = save_dir + "long_binary_evolution_plot_" + str(binary_tags[0]) + "_" + str(binary_tags[1])
                
                plt.clf()
                fig = plt.figure()
                fig.set_size_inches(6, 14.)
                gs = gridspec.GridSpec(6, 1)
                gs.update(hspace=0.0)
                ax1 = fig.add_subplot(gs[0,0]) #Separation
                ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
                ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #Eccentricity
                ax4 = fig.add_subplot(gs[3,0], sharex=ax1)
                ax5 = fig.add_subplot(gs[4,0], sharex=ax1)
                ax6 = fig.add_subplot(gs[5,0], sharex=ax1)

                ax1.semilogy(particle_data['time'].value, separation.value, lw=0.5)#, color=line_colour[nit])
                if args.plot_matched_times != 'False':
                    matches_dict = {}
                    with open('./quadrupole_rv_times.csv','rU') as match_times:
                        for line in match_times:
                            row = re.split('[[]|[]]', line)
                            if row[0][0] != '#':
                                no_thres = eval('['+row[1]+']')
                                obs_thres = eval('['+row[3]+']')
                                matches_dict.update({str(row[0][0]):[no_thres, obs_thres]})
                    proj_colours = plt.cm.nipy_spectral(np.linspace(0,1,8))
                    #proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']
                    for proj in range(8):
                        alpha = [0.5, 1.0]
                        for thres_it in range(2):
                            for match_times in matches_dict[str(proj)][thres_it]:
                                try:
                                    ax1.axvspan(match_times[0], match_times[1], alpha=thres_it, color=proj_colours[proj])
                                except:
                                    ax1.axvline(match_times, alpha=thres_it, color=proj_colours[proj], lw=2)
                    
                ax1.set_ylabel('Separation (AU)')
                ax1.axhline(y=25, color='k')
                plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
                if args.end_time != None:
                    ax1.set_xlim(right=float(args.end_time))
                else:
                    ax1.set_xlim([particle_data['time'][non_nan_inds][0].value, particle_data['time'][non_nan_inds][-1].value])
                ax1.yaxis.set_ticks_position('both')
                ax1.tick_params(axis='y', which='major', direction="in")
                ax1.tick_params(axis='y', which='minor', direction="in")
                ax1.tick_params(axis='x', which='major', direction="in")

                #Plot accretion
                for pit in range(len(labels)):
                    ax2.semilogy(particle_data['time'].value, binary_accretion[pit].value, label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
                total_accretion = np.sum(binary_accretion, axis=0)
                ax2.legend(loc='best')
                ax2.semilogy(particle_data['time'].value, total_accretion.value, color='k', label="Total", lw=0.5)
                ax2.set_ylabel(r'Accretion Rate (M$_\odot$/yr)')
                ax2.yaxis.set_ticks_position('both')
                ax2.tick_params(axis='y', which='major', direction="in")
                ax2.tick_params(axis='y', which='minor', direction="in")
                ax2.tick_params(axis='x', which='major', direction="in")
                plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)

                for pit in range(len(labels)):
                    ax3.semilogy(particle_data['time'].value, binary_masses[pit].value, label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
                total_mass = np.sum(binary_masses, axis=0)
                ax3.semilogy(particle_data['time'].value, total_mass.value, color='k', label="Total", lw=0.5)
                ax3.set_ylabel(r'Mass (M$_\odot$)')
                ax3.set_ylim(bottom=1.e-2)
                ax3.yaxis.set_ticks_position('both')
                ax3.tick_params(axis='y', which='major', direction="in")
                ax3.tick_params(axis='y', which='minor', direction="in")
                ax3.tick_params(axis='x', which='major', direction="in")
                plt.setp([ax3.get_xticklabels() for ax3 in fig.axes[:-1]], visible=False)
                
                ax4.semilogy(particle_data['time'].value, e.value, lw=0.5)
                ax4.set_ylabel(r'Eccentricity')
                ax4.yaxis.set_ticks_position('both')
                ax4.tick_params(axis='y', which='major', direction="in")
                ax4.tick_params(axis='y', which='minor', direction="in")
                ax4.tick_params(axis='x', which='major', direction="in")
                
                ax5.semilogy(particle_data['time'].value, (abs(E_pot)/total_mass).value, lw=0.5, label='Potential')
                ax5.semilogy(particle_data['time'].value, (E_kin.value/total_mass), lw=0.5, label='Kinetic')
                ax5.legend(loc='best')
                ax5.set_ylabel(r'Spec. En. (erg/g)')
                ax5.yaxis.set_ticks_position('both')
                ax5.tick_params(axis='y', which='major', direction="in")
                ax5.tick_params(axis='y', which='minor', direction="in")
                ax5.tick_params(axis='x', which='major', direction="in")
                
                ax6.semilogy(particle_data['time'].value, (np.sqrt(np.sum(L_1**2, axis=0))/total_mass).value, alpha=0.5, lw=0.5, label=labels[0])
                ax6.semilogy(particle_data['time'].value, (np.sqrt(np.sum(L_2**2, axis=0))/total_mass).value, alpha=0.5, lw=0.5, label=labels[1])
                ax6.semilogy(particle_data['time'].value, (L_tot/total_mass).value, lw=0.5, label='Total')
                ax6.legend(loc='best')
                ax6.set_ylabel(r'Spec. Ang. Mom. ($cm^2/s$)')
                ax6.yaxis.set_ticks_position('both')
                ax6.tick_params(axis='y', which='major', direction="in")
                ax6.tick_params(axis='y', which='minor', direction="in")
                ax6.tick_params(axis='x', which='major', direction="in")
                ax6.set_xlabel('Time (yr)')
                
                
                #Save image
                plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
                print("Created image", image_name)
                
    
#Save reduced system data:
file = open(save_dir+'reduced_system_data.pkl', 'wb')
pickle.dump((reduced_systems_data), file)
file.close()
