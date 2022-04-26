import numpy as np
import yt
import pickle
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import collections

#Define globals
f_acc= 0.5
Accretion_array = []

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-verbose", "--verbose_printing", help="Would you like to print debug lines?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    

#=====================================================================================================

args = parse_inputs()

#==========================================================================================

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]

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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
    
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()

try:
    file = open("sink_birth_conditions.pkl", 'wb')
    Sink_bound_birth = pickle.load(file)
    file.close()
except:
    Sink_bound_birth = []
Mass_plus_blank_row = np.vstack([np.zeros(len(global_data['m'][0])), global_data['m']])
diff_arr =  (Mass_plus_blank_row[1:]-Mass_plus_blank_row[:-1])
zero_inds = np.where(diff_arr == 0)
diff_arr[zero_inds] = 1
formation_inds = np.where(diff_arr == global_data['m'])
if len(Sink_bound_birth) > 0:
    import pdb
    pdb.set_trace()
n_stars = 0
for sink_id in formation_inds[1]:
    form_time_it = formation_inds[0][sink_id]
    new_sink_pos = np.array([global_data['x'][form_time_it][sink_id], global_data['y'][form_time_it][sink_id], global_data['z'][form_time_it][sink_id]]).T
    new_sink_vel = np.array([global_data['ux'][form_time_it][sink_id], global_data['uy'][form_time_it][sink_id], global_data['uz'][form_time_it][sink_id]]).T
    new_sink_mass = np.array(global_data['m'][form_time_it][sink_id])

    abspos = np.array([global_data['x'][form_time_it][:sink_id], global_data['y'][form_time_it][:sink_id], global_data['z'][form_time_it][:sink_id]]).T
    absvel = np.array([global_data['ux'][form_time_it][:sink_id], global_data['uy'][form_time_it][:sink_id], global_data['uz'][form_time_it][:sink_id]]).T
    mass = np.array(global_data['m'][form_time_it][:sink_id])
    
    rel_pos = abspos - new_sink_pos
    update_seps_neg = np.argwhere(abs(rel_pos)<-0.5)
    update_seps_pos = np.argwhere(abs(rel_pos)>0.5)
    import pdb
    pdb.set_trace()
    rel_pos[update_seps_neg] = rel_pos[update_seps_neg] + 0.5
    rel_pos[update_seps_pos] = rel_pos[update_seps_pos] - 0.5
    rel_sep = np.sqrt(rel_pos[:,0]**2 + rel_pos[:,1]**2 + rel_pos[:,2]**2)
    rel_vel = absvel - new_sink_vel
    rel_speed = np.sqrt(rel_vel[:,0]**2 + rel_vel[:,1]**2 + rel_vel[:,2]**2)
    mtm = new_sink_mass * mass
    mpm = new_sink_mass + mass
    newtonianPotential = -1./rel_sep
    
    Ekin = 0.5 * mtm/mpm * rel_speed**2
    Epot = Grho * mtm * newtonianPotential
    Etot = Ekin + Epot
    
    if True in (Etot<0):
        born_bound = True
        most_bound_sink_id = np.argmin(Etot)
        lowest_Etot = np.nanmin(Etot)
        #Do multiplicity analysis
        time_it = formation_inds[0][sink_id]
        n_stars = np.where(global_data['m'][time_it]>0)[0]

        abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
        absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
        mass = np.array(global_data['m'][time_it][n_stars])
        time = global_data['time'][time_it][n_stars][0]
        S = pr.Sink()
        S._jet_factor = 1.
        S._scale_l = scale_l.value
        S._scale_v = scale_v.value
        S._scale_t = scale_t.value
        S._scale_d = scale_d.value
        S._time = yt.YTArray(time, '')
        S._abspos = yt.YTArray(abspos, '')
        S._absvel = yt.YTArray(absvel, '')
        S._mass = yt.YTArray(mass, '')
        res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho)
        multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
        most_bound_sep = np.nan
        for multi_ind in multi_inds:
            if np.isnan(most_bound_sep):
                sys_comps = losi(multi_ind, res)
                sys_string = sorted(flatten(sys_comps))
                if sink_id in sys_string:
                    sys_comps = str(sys_comps)
                    if sink_id in sys_string:
                        reduced = False
                        while reduced == False:
                            open_braket_ind = []
                            for char_it in range(len(sys_comps)):
                                if sys_comps[char_it] == '[':
                                    open_braket_ind.append(char_it)
                                if sys_comps[char_it] == ']':
                                    open_ind = open_braket_ind.pop()
                                    sub_sys = eval(sys_comps[open_ind:char_it+1])
                                    if len(sub_sys) == 2:
                                        ind_1 = sub_sys[0]
                                        ind_2 = sub_sys[1]
                                        pos_diff = res['abspos'][ind_1] - res['abspos'][ind_2]
                                        sep_value = np.sqrt(np.sum(pos_diff**2))
                                        position = yt.YTArray(res['abspos'][np.array([ind_1,ind_2])].T, 'au')
                                        if sep_value > (scale_l.in_units('AU')/2):
                                            update_inds = np.where(abs(pos_diff)>scale_l.in_units('AU')/2)[0]
                                            for ind in update_inds:
                                                if pos_diff[ind] < 0:
                                                    pos_diff[ind] = pos_diff[ind] + scale_l.in_units('AU').value
                                                else:
                                                    pos_diff[ind] = pos_diff[ind] - scale_l.in_units('AU').value
                                            sep_value = np.sqrt(np.sum(pos_diff**2))
                                            for ind in update_inds:
                                                pos_update_ind = np.where(position[ind]<scale_l.in_units('AU')/2)[0]
                                                position[ind][pos_update_ind] = position[ind][pos_update_ind] + scale_l.in_units('AU')
                                                #velocity[ind][pos_update_ind] = -1*velocity[ind][pos_update_ind]# + scale_l.in_units('AU')
                                            pos_diff = position[ind][0] - position[ind][0]
                                            sep_value = np.sqrt(np.sum(pos_diff**2))
                                            if sep_value > 20000:
                                                import pdb
                                                pdb.set_trace()
                                        if sink_id in sub_sys:
                                            reduced = True
                                            most_bound_sep = sep_value
                                            first_bound_sink = sub_sys[np.argwhere(sub_sys != sink_id)[0][0]]
                                            first_bound_sink = losi(first_bound_sink, res)
                                            if most_bound_sink_id != first_bound_sink:
                                                import pdb
                                                pdb.set_trace()
                                            break
                                        replace_ind = np.where((res['index1']==sub_sys[0])&(res['index2']==sub_sys[1]))[0][0]
                                        replace_string = str(replace_ind)
                                        sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                                        if '[' not in sys_comps:
                                            reduced = True
                                        elif '[' in sys_comps:
                                            if len(eval(sys_comps)) == 1:
                                                reduced = True
                                        break
                                    else:
                                        import pdb
                                        pdb.set_trace()
    else:
        born_bound = False
        if len(Etot) > 0:
            most_bound_sink_id = np.argmin(Etot)
        else:
            most_bound_sink_id = np.nan
        most_bound_sep = np.nan
        first_bound_sink = np.nan
        lowest_Etot = np.nan
        
        time_it = formation_inds[0][sink_id]
        new_sink_pos_x = global_data['x'][time_it:,sink_id]
        new_sink_pos_y = global_data['y'][time_it:,sink_id]
        new_sink_pos_z = global_data['z'][time_it:,sink_id]
        
        new_sink_vel_x = global_data['ux'][time_it:,sink_id]
        new_sink_vel_y = global_data['uy'][time_it:,sink_id]
        new_sink_vel_z = global_data['uz'][time_it:,sink_id]
        
        new_sink_mass = np.array(global_data['m'][time_it:,sink_id])
        
        abspos_x = global_data['x'][time_it:]
        abspos_y = global_data['y'][time_it:]
        abspos_z = global_data['z'][time_it:]
        
        absvel_x = global_data['ux'][time_it:]
        absvel_y = global_data['uy'][time_it:]
        absvel_z = global_data['uz'][time_it:]
        
        mass = np.array(global_data['m'][time_it:])
        
        rel_pos_x = abspos_x.T - new_sink_pos_x
        rel_pos_y = abspos_y.T - new_sink_pos_y
        rel_pos_z = abspos_z.T - new_sink_pos_z
        import pdb
        pdb.set_trace()
        update_seps_x_neg = np.argwhere(abs(rel_pos_x)>-0.5)
        update_seps_x_pos = np.argwhere(abs(rel_pos_x)<0.5)
        update_seps = np.array(np.where(abs(rel_pos)>0.5)).T
        for update_sep in update_seps:
            if rel_pos[update_sep[0]][update_sep[1]][update_sep[2]] < 0:
                rel_pos[update_sep[0]][update_sep[1]][update_sep[2]] = rel_pos[update_sep[0]][update_sep[1]][update_sep[2]] + 0.5
            else:
                rel_pos[update_sep[0]][update_sep[1]][update_sep[2]] = rel_pos[update_sep[0]][update_sep[1]][update_sep[2]] - 0.5
        rel_sep = np.sqrt(np.sum(np.square(rel_pos), axis=2))
        rel_vel = absvel - new_sink_vel
        rel_speed = np.sqrt(np.sum(np.square(rel_vel), axis=2))
        mtm = new_sink_mass * mass.T
        mpm = new_sink_mass + mass.T
        newtonianPotential = -1./rel_sep
        Ekin = 0.5 * mtm/mpm * rel_speed**2
        Epot = Grho * mtm * newtonianPotential
        Etot = Ekin + Epot
        Etot[Etot == -1*np.inf] = np.nan
        Etot_min = np.nanmin(Etot, axis=0)
        if True in (Etot_min<0):
            time_it = formation_inds[0][sink_id] + np.where(Etot_min<0)[0][0]
            lowest_Etot = Etot_min[np.where(Etot_min<0)[0][0]]
            n_stars = np.where(global_data['m'][time_it]>0)[0]
            
            abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
            absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
            mass = np.array(global_data['m'][time_it][n_stars])
            time = global_data['time'][time_it][n_stars][0]
            S = pr.Sink()
            S._jet_factor = 1.
            S._scale_l = scale_l.value
            S._scale_v = scale_v.value
            S._scale_t = scale_t.value
            S._scale_d = scale_d.value
            S._time = yt.YTArray(time, '')
            S._abspos = yt.YTArray(abspos, '')
            S._absvel = yt.YTArray(absvel, '')
            S._mass = yt.YTArray(mass, '')
            res = m.multipleAnalysis(S,cutoff=10000, bound_check=True, nmax=6, cyclic=True, Grho=Grho)
            multi_inds = np.where((res['n']>1) & (res['topSystem']==True))[0]
            most_bound_sep = np.nan
            for multi_ind in multi_inds:
                if np.isnan(most_bound_sep):
                    sys_comps = losi(multi_ind, res)
                    sys_string = sorted(flatten(sys_comps))
                    if sink_id in sys_string:
                        sys_comps = str(sys_comps)
                        if sink_id in sys_string:
                            reduced = False
                            while reduced == False:
                                open_braket_ind = []
                                for char_it in range(len(sys_comps)):
                                    if sys_comps[char_it] == '[':
                                        open_braket_ind.append(char_it)
                                    if sys_comps[char_it] == ']':
                                        open_ind = open_braket_ind.pop()
                                        sub_sys = eval(sys_comps[open_ind:char_it+1])
                                        if len(sub_sys) == 2:
                                            ind_1 = sub_sys[0]
                                            ind_2 = sub_sys[1]
                                            pos_diff = res['abspos'][ind_1] - res['abspos'][ind_2]
                                            sep_value = np.sqrt(np.sum(pos_diff**2))
                                            position = yt.YTArray(res['abspos'][np.array([ind_1,ind_2])].T, 'au')
                                            if sep_value > (scale_l.in_units('AU')/2):
                                                update_inds = np.where(abs(pos_diff)>scale_l.in_units('AU')/2)[0]
                                                for ind in update_inds:
                                                    if pos_diff[ind] < 0:
                                                        pos_diff[ind] = pos_diff[ind] + scale_l.in_units('AU').value
                                                    else:
                                                        pos_diff[ind] = pos_diff[ind] - scale_l.in_units('AU').value
                                                sep_value = np.sqrt(np.sum(pos_diff**2))
                                                for ind in update_inds:
                                                    pos_update_ind = np.where(position[ind]<scale_l.in_units('AU')/2)[0]
                                                    position[ind][pos_update_ind] = position[ind][pos_update_ind] + scale_l.in_units('AU')
                                                    #velocity[ind][pos_update_ind] = -1*velocity[ind][pos_update_ind]# + scale_l.in_units('AU')
                                                pos_diff = position[ind][0] - position[ind][0]
                                                sep_value = np.sqrt(np.sum(pos_diff**2))
                                                if sep_value > 20000:
                                                    import pdb
                                                    pdb.set_trace()
                                            if sink_id in sub_sys:
                                                reduced = True
                                                most_bound_sep = sep_value
                                                first_bound_sink = sub_sys[np.argwhere(sub_sys != sink_id)[0][0]]
                                                first_bound_sink = losi(first_bound_sink, res)
                                                break
                                            replace_ind = np.where((res['index1']==sub_sys[0])&(res['index2']==sub_sys[1]))[0][0]
                                            replace_string = str(replace_ind)
                                            sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                                            if '[' not in sys_comps:
                                                reduced = True
                                            elif '[' in sys_comps:
                                                if len(eval(sys_comps)) == 1:
                                                    reduced = True
                                            break
                                        else:
                                            import pdb
                                            pdb.set_trace()
    if np.isnan(most_bound_sep) and lowest_Etot < 0:
        import pdb
        pdb.set_trace()
    Sink_bound_birth.append([born_bound, most_bound_sink_id, first_bound_sink, most_bound_sep, lowest_Etot])
    print("Birth conditions of sink", sink_id, "is", Sink_bound_birth[-1])

    file = open("sink_birth_conditions.pkl", 'wb')
    pickle.dump((Sink_bound_birth),file)
    file.close()
