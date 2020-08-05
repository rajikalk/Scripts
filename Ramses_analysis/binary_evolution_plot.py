import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pyramses import rsink
import sys
import os
import yt
from scipy.stats import mode
import numpy.ma as ma

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-radius", "--proximity_radius", help="within what radius (in AU) do you want to save sink particle data?", type=float, default=10000.0)
    parser.add_argument("-def_sys", "--define_system", help="Is there a particular system that you would like to plot?", type=str, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def boolean_indexing(v, fillval=np.nan):
    """
    function for squaring off array, so that they are retangular arrays. It filled empty indexes with NaNs by default.
    """
    lens = np.array([len(item) for item in v])
    mask = lens[:,None] > np.arange(lens.max())
    out = np.full(mask.shape,fillval)
    out[mask] = np.concatenate(v)
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
    
units = {"length_unit":yt.YTQuantity(4.0,"pc"), "mass_unit":yt.YTQuantity(2998,"Msun"), "velocity_unit":yt.YTQuantity(0.18, "km/s"), "time_unit":yt.YTQuantity(685706129102738.9, "s"), "density_unit":yt.YTQuantity(46.84375, "Msun/pc**3")}

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
            loaded_sink_data = loaded_sink_data[counter:]
            if counter < len(loaded_sink_data):
                updating = True
                print('pickle data is not up to date! Updating')
        except:
            os.system('cp '+save_dir+'particle_data_raw_tmp.pkl '+save_dir+'particle_data_raw.pkl ')
            file_open = open(save_dir+'particle_data_raw.pkl', 'rb')
            particle_data, counter, sink_ind = pickle.load(file_open)
            file_open.close()
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
        particle_data.update({'separation':[]})
        particle_data.update({'eccentricity':[]})
        particle_data.update({'L_orb':[]})
        particle_data.update({'L_orb_tot':[]})
        particle_data.update({'dist_from_CoM':[]})
        counter = 0
        
    sink_form_time = 0

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
                velocity = yt.YTArray([particle_data['velx'][-1],particle_data['vely'][-1], particle_data['velz'][-1]], 'au')
            
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

                    separation = np.sum(distance_from_com)
                    relative_speed_to_com = np.sqrt(np.sum(vel_rel_to_com**2, axis=0))
                
                    #Calculmate mu and orbital energy
                    reduced_mass = np.product(particle_data['mass'][-1].in_units('g'))/np.sum(particle_data['mass'][-1].in_units('g'))
                    E_pot = (-1*(yt.units.G*np.product(particle_data['mass'][-1].in_units('g')))/separation.in_units('cm')).in_units('erg')
                    E_kin = np.sum((0.5*particle_data['mass'][-1].in_units('g')*relative_speed_to_com.in_units('cm/s')**2).in_units('erg'))
                    
                    #epsilon is the specific orbital energy
                    epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                    
                    #Calculate orbital energy
                    m_x_r = yt.YTArray(np.cross(pos_rel_to_com.T.in_units('cm'), vel_rel_to_com.T.in_units('cm/s')).T, 'cm**2/s')
                    L = particle_data['mass'][-1].in_units('g').T*m_x_r
                    L_tot = np.sqrt(np.sum(np.sum(L, axis=1)**2, axis=0))
                    
                    #h_val is the specific angular momentum
                    h_val = L_tot/reduced_mass.in_units('g')
                    
                    e = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.G*np.sum(particle_data['mass'][-1].in_units('g')))**2.))
                    particle_data['eccentricity'].append(e)
                    particle_data['separation'].append(separation.in_units('au'))
                    particle_data['L_orb'].append(L)
                    particle_data['L_orb_tot'].append(L_tot)
                except:
                    particle_data['eccentricity'].append(yt.YTArray(np.nan, ''))
                    particle_data['separation'].append(yt.YTArray(np.nan, 'au'))
                    particle_data['L_orb'].append(yt.YTArray([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]], 'cm**2*g/s'))
                    particle_data['L_orb_tot'].append(yt.YTArray([np.nan, np.nan, np.nan], 'cm**2*g/s'))
                """
                if len(non_nan_inds) == 3:
                    #Calculate how strongly bound the other particles are to the primary
                    #get positions
                    position = yt.YTArray([particle_data['posx'][-1][non_nan_inds],particle_data['posy'][-1][non_nan_inds], particle_data['posz'][-1][non_nan_inds]], 'au')
                    velocity = yt.YTArray([particle_data['velx'][-1][non_nan_inds],particle_data['vely'][-1][non_nan_inds], particle_data['velz'][-1][non_nan_inds]], 'au')
                    
                    temporary_potentials = []
                    most_bound_partical_ind = []
                    for particle_it in range(len(non_nan_inds)):
                        current_pos = position[particle_it]
                        current_vel = velocity[particle_it]
                        current_mass = particle_data['mass'][-1][particle_it]
                        separations = np.sqrt(np.sum((current_pos.T - position)**2, axis=1))
                        relative_speed = np.sqrt(np.sum((current_vel.T - velocity)**2, axis=1))
                        
                        E_pot = (-1*(yt.units.G*current_mass.in_units('g')*particle_data['mass'][-1][non_nan_inds].in_units('g'))/separations.in_units('cm')).in_units('erg')
                        #E_kin = (0.5*particle_data['mass'][-1][non_nan_inds].in_units('g')*relative_speed.in_units('cm/s')**2).in_units('erg')
                        #boundness = E_pot + E_kin
                        most_bound_part_ind = np.argsort(E_pot)[1]
                        most_bound_partical_ind.append(particle_data['particle_tag'][most_bound_part_ind])
                        temporary_potentials.append(E_pot[most_bound_part_ind])
                        
                    binary_potential = mode(temporary_potentials)[0][0]
                    hierarchical_inds = np.argwhere(temporary_potentials == binary_potential).T[0].tolist()
                    for ind in range(len(non_nan_inds)):
                        if ind not in hierarchical_inds:
                            hierarchical_inds.append(ind)
                    
                    #now iterate to calculate central binary properties, and then update to calculate trinary properties
                    position = [[particle_data['posx'][-1][hierarchical_inds[0]]],[particle_data['posy'][-1][hierarchical_inds[0]]], [particle_data['posz'][-1][hierarchical_inds[0]]]]
                    velocity = [[particle_data['velx'][-1][hierarchical_inds[0]]],[particle_data['vely'][-1][hierarchical_inds[0]]], [particle_data['velz'][-1][hierarchical_inds[0]]]]
                    mass = [particle_data['mass'][-1][hierarchical_inds[0]]]
                    for hit in range(1,2):
                        position[0].append(particle_data['posx'][-1][hierarchical_inds[hit]])
                        position[1].append(particle_data['posy'][-1][hierarchical_inds[hit]])
                        position[2].append(particle_data['posz'][-1][hierarchical_inds[hit]])
                        
                        velocity[0].append(particle_data['velx'][-1][hierarchical_inds[hit]])
                        velocity[1].append(particle_data['vely'][-1][hierarchical_inds[hit]])
                        velocity[2].append(particle_data['velz'][-1][hierarchical_inds[hit]])
                        
                        mass.append(particle_data['mass'][-1][hierarchical_inds[hit]])
                        import pdb
                        pdb.set_trace()
                        
                        CoM_pos
                    
                    if len(non_nan_inds) > 2:
                        import pdb
                        pdb.set_trace()
                    '''
                    position = yt.YTArray([particle_data['posx'][-1],particle_data['posy'][-1], particle_data['posz'][-1]], 'au')
                    velocity = yt.YTArray([particle_data['velx'][-1],particle_data['vely'][-1], particle_data['velz'][-1]], 'au')
                   
                    CoM_pos = np.sum((position.T[non_nan_inds].T*particle_data['mass'][-1][non_nan_inds]).T, axis=0)/np.sum(particle_data['mass'][-1][non_nan_inds])
                    CoM_vel = np.sum((velocity.T[non_nan_inds].T*particle_data['mass'][-1][non_nan_inds]).T, axis=0)/np.sum(particle_data['mass'][-1][non_nan_inds])
                   
                    pos_rel_to_com = (position.T - CoM_pos).T
                    vel_rel_to_com = (velocity.T - CoM_vel).T
                    distance_from_com = np.sqrt(np.sum(pos_rel_to_com**2, axis=0))
                    relative_speed = np.sqrt(np.sum(vel_rel_to_com**2, axis=0))
                    particle_data['separation'].append([distance_from_com.in_units('au')])
                   
                    #Calculmate mu
                    mass_total = np.sum(particle_data['mass'][-1][non_nan_inds])
                    #mu = yt.units.G*mass_total.in_units('g')*particle_data['mass'][-1].in_units('g')
                    mu = yt.units.G*mass_total.in_units('g')
                    epsilon = (relative_speed.in_units('cm/s')**2)/2 - mu/distance_from_com.in_units('cm')
                    
                    reduced_mass = yt.YTArray(np.product(particle_data['mass'][-1][non_nan_inds])/np.sum(particle_data['mass'][-1][non_nan_inds]), 'msun')
                    m_x_r = yt.YTArray(np.cross(pos_rel_to_com.T.in_units('cm'), vel_rel_to_com.T.in_units('cm/s')).T, 'cm**2/s')
                    L = particle_data['mass'][-1].in_units('g').T*m_x_r
                    particle_data['L_orb'].append([L])
                    L_tot = np.sum(L**2, axis=0)
                    particle_data['L_orb_tot'].append([L_tot])
                    h_val = np.sqrt(np.sum(L_tot[non_nan_inds]))/reduced_mass.in_units('g')
                    
                    e = np.sqrt(1 + (2.*epsilon*h_val**2.)/(mu**2.))
                    particle_data['eccentricity'].append([e])
                    '''
                else:
                    print("Not a binary!")
                """

    print("Finished reading particle data")

print('Reading pickle file')
file_open = open(save_dir+'particle_data_raw.pkl', 'rb')
particle_data, counter, sink_ind = pickle.load(file_open)
file_open.close()

print("making arrays neat")

import pdb
pdb.set_trace()

for key in list(particle_data.keys()):
    try:
        particle_data[key] = boolean_indexing(particle_data[key])
    except:
        print("1D array can't be squared off")

try:
    for key in list(particle_data.keys()):
        print("key =", key)
        try:
            unit_string = particle_data[key][0].units
            print('unit_string =', unit_string)
        except:
            unit_string = ""
            print('unit_string =', unit_string)
            
        particle_data[key] = yt.YTArray(np.array(particle_data[key]), unit_string)
except:
    print('arrays have already been converted to the YT Array')


'''
image_name = save_dir + "separation_resolution_study"
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
    ax1.semilogy(particle_data['time'].value, particle_data['dist_from_CoM'].T[nit].value, label=labels[nit], lw=0.5)#, color=line_colour[nit])
ax1.set_ylabel('Distance from CoM (AU)')
ax1.legend(loc='upper left')
ax1.set_ylim(bottom=1.e0)
plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
ax1.set_xlim([particle_data['time'][0].value, particle_data['time'][-1].value])

#Plot accretion
for pit in range(len(labels)):
    ax2.semilogy(particle_data['time'].value, particle_data['mdot'].value.T[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
total_accretion = np.sum(np.nan_to_num(particle_data['mdot'].value.T), axis=0)
ax2.semilogy(particle_data['time'].value, total_accretion, color='k', label="Total", lw=0.5)
ax2.set_ylabel(r'Accretion Rate (M$_\odot$/yr)')
#ax2.legend(loc='best')
plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)

for pit in range(len(labels)):
    ax3.semilogy(particle_data['time'].value, particle_data['mass'].value.T[pit], label=labels[pit], alpha=0.5, lw=0.5)#, color=line_colour[pit])
total_accretion = np.sum(np.nan_to_num(particle_data['mass'].value.T), axis=0)
ax3.semilogy(particle_data['time'].value, total_accretion, color='k', label="Total", lw=0.5)
ax3.set_ylabel(r'Mass (M$_\odot$)')
ax3.set_ylim(bottom=1.e-2)
#ax2.legend(loc='best')

#Save image
plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
print("Created image", image_name)
'''

#Anaylsys particular systems
if args.define_system != None:
    systems_hierarchy = args.define_system
else:
    print("Not making plots for the multiple star system")

if systems_hierarchy != None:
    while len(find(systems_hierarchy, '[')) > 1:
        indexes = []
        for ind, char in enumerate(systems_hierarchy):
            if char == '[':
                indexes.append(ind)
            if char == ']':
                start_ind = indexes.pop()
                end_ind = ind
                
                binary_tags = eval('['+systems_hierarchy[start_ind+1:end_ind]+']')
                system_tag = int(np.mean(binary_tags))
                
                #Found binary, now to create  necessary plots
                particle_inds = []
                for tag in binary_tags:
                    try:
                        ind = np.argwhere(particle_data['particle_tag'] == tag)[0][0]
                    except:
                        ind = particle_data['particle_tag'].index(tag)
                    particle_inds.append(ind)
                    
                import pdb
                pdb.set_trace()
                    
                #define array like position, velocity, mass
                binary_masses = particle_data['mass'].T[particle_inds]
                binary_positions = np.array([particle_data['posx'].T[particle_inds],particle_data['posy'].T[particle_inds], particle_data['posz'].T[particle_inds]])
                binary_velocities = np.array([particle_data['velx'].T[particle_inds],particle_data['vely'].T[particle_inds], particle_data['velz'].T[particle_inds]])
                
                CoM_pos = np.sum((binary_positions*binary_masses).T, axis=0)/np.sum(binary_masses)
                CoM_vel = np.sum((binary_velocities*binary_masses).T, axis=0)/np.sum(binary_masses)
                
                pos_rel_to_com = (position.T - CoM_pos).T
                distance_from_com = np.sqrt(np.sum(pos_rel_to_com**2, axis=0))
                
                vel_rel_to_com = (velocity.T - CoM_vel).T

                separation = np.sum(distance_from_com)
                relative_speed_to_com = np.sqrt(np.sum(vel_rel_to_com**2, axis=0))
            
                #Calculmate mu and orbital energy
                reduced_mass = np.product(binary_masses.in_units('g'))/np.sum(binary_masses.in_units('g'))
                E_pot = (-1*(yt.units.G*np.product(binary_masses.in_units('g')))/separation.in_units('cm')).in_units('erg')
                E_kin = np.sum((0.5*binary_masses.in_units('g')*relative_speed_to_com.in_units('cm/s')**2).in_units('erg'))
                
                #epsilon is the specific orbital energy
                epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
                
                #Calculate orbital energy
                m_x_r = yt.YTArray(np.cross(pos_rel_to_com.T.in_units('cm'), vel_rel_to_com.T.in_units('cm/s')).T, 'cm**2/s')
                L = particle_data['mass'][-1].in_units('g').T*m_x_r
                L_tot = np.sqrt(np.sum(np.sum(L, axis=1)**2, axis=0))
                
                #h_val is the specific angular momentum
                h_val = L_tot/reduced_mass.in_units('g')
                
                e = np.sqrt(1 + (2.*epsilon*h_val**2.)/((yt.units.G*np.sum(particle_data['mass'][-1].in_units('g')))**2.))
                
                systems_hierarchy = systems_hierarchy[:start_ind] + str(system_tag) + systems_hierarchy[end_ind+1:]
                print('systems_hierarchy =', systems_hierarchy)
                break
    
    

