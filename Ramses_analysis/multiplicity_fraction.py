import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
import sys
import collections

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-proj", "--projected_separation", help="do you want to use projected separation instead of true separation?", default="False", type=str)
    parser.add_argument("-ax", "--axis", help="what axis do you want to project separation onto?", default='x', type=str)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
def losi(i, res):
    if (res['n'][i]==1):
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
        
#=====================================================================================================

#Baraffe models
        
#=====================================================================================================

args = parse_inputs()

units = {"length_unit":yt.YTQuantity(4.0,"pc"), "mass_unit":yt.YTQuantity(2998,"Msun"), "velocity_unit":yt.YTQuantity(0.18, "km/s"), "time_unit":yt.YTQuantity(685706129102738.9, "s"), "density_unit":yt.YTQuantity(46.84375, "Msun/pc**3")}

scale_l = 1.23427103e19 # 4 pc
scale_v = 1.8e4         # 0.18 km/s == sound speed
scale_t = 6.85706128e14 # 4 pc / 0.18 km/s
scale_d = 3.171441e-21  # 2998 Msun / (4 pc)^3

'''
#low resolution data
datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data_256'
nout = 133
'''

#high resolution data
#datadir = '/lustre/astro/troels/IMF_512/binary_analysis/data'
datadir = sys.argv[1]
savedir = sys.argv[2]
file_no_range = [298,407]


Separations = []
Times = []
Multiplicities = []
Companion_frequencies = []
N_sys_total = []
CF_Array_Full = []
prev_time = np.nan
prev_inds = []
prev_components = []
prev_mass = []
curr_time = np.nan
curr_inds = []
curr_components = []
curr_mass = []
curr_res = np.nan
next_time = np.nan
next_inds = []
next_components = []
next_mass = []

for nout in range(file_no_range[0]-1, file_no_range[1]+1):
    # load sink data from snapshot nout in to S
    S = pr.Sink()
    S.nout = nout
    S.datadir = datadir
    try:
        #Get data from previous time step and next time step:
        S.rsink()
        S._jet_factor = 1.
        S._scale_l = scale_l
        S._scale_v = scale_v
        S._scale_t = scale_t
        S._scale_d = scale_d
        i = S.read_info_file(nout,datadir=datadir)
        S._time = i['time']
        res = m.multipleAnalysis(S,nmax=len(S.mass),cutoff=1e4)
        if np.isnan(prev_time):
            # do multiplicity analysis
            #Cutoff of 1e4 AU is because this is the cutoff that Tobin et al (2016) uses.

            # find all single systems and all multiple systems
            s = np.where(res['n']==1)[0]
            b = np.where(res['n']==2)[0]
            t = np.where(res['n']==3)[0]
            q = np.where(res['n']==4)[0]
            q5 = np.where(res['n']==5)[0]
            s6 = np.where(res['n']==6)[0]
            h = np.where(res['n']==7)[0]
            o = np.where(res['n']==8)[0]
            multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist() + h.tolist() + o.tolist())
            
            #save components
            for multi_ind in multi_inds:
                sys_comps = losi(multi_ind, res)
                sys_comps = sorted(flatten(sys_comps))
                prev_components.append(sys_comps)
                prev_mass.append(res['mass'][multi_ind])
                prev_inds.append(multi_ind)
            prev_time = S._time * units['time_unit'].in_units('yr')
            print('set initial prev_time')
        elif np.isnan(curr_time):
            curr_res = res
            
            s = np.where(res['n']==1)[0]
            b = np.where(res['n']==2)[0]
            t = np.where(res['n']==3)[0]
            q = np.where(res['n']==4)[0]
            q5 = np.where(res['n']==5)[0]
            s6 = np.where(res['n']==6)[0]
            h = np.where(res['n']==7)[0]
            o = np.where(res['n']==8)[0]
            multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist() + h.tolist() + o.tolist())
            
            #save components
            for multi_ind in multi_inds:
                sys_comps = losi(multi_ind, res)
                sys_comps = sorted(flatten(sys_comps))
                curr_components.append(sys_comps)
                curr_mass.append(res['mass'][multi_ind])
                curr_inds.append(multi_ind)
            curr_time = S._time * units['time_unit'].in_units('yr')
            print('set initial curr_time')
        else:
            next_time = np.nan
            next_inds = []
            next_components = []
            next_mass = []
            
            next_time = S._time * units['time_unit'].in_units('yr')
            
            s = np.where(res['n']==1)[0]
            b = np.where(res['n']==2)[0]
            t = np.where(res['n']==3)[0]
            q = np.where(res['n']==4)[0]
            q5 = np.where(res['n']==5)[0]
            s6 = np.where(res['n']==6)[0]
            h = np.where(res['n']==7)[0]
            o = np.where(res['n']==8)[0]
            multi_inds = np.array(b.tolist() + t.tolist() + q.tolist() + q5.tolist() + s6.tolist() + h.tolist() + o.tolist())
            
            #save components
            for multi_ind in multi_inds:
                sys_comps = losi(multi_ind, res)
                sys_comps = sorted(flatten(sys_comps))
                next_components.append(sys_comps)
                next_mass.append(res['mass'][multi_ind])
                next_inds.append(multi_ind)
                    
            #CALCULATE MAGNITUDE
            radius = yt.YTQuantity(2.0, 'rsun')
            temperature = yt.YTQuantity(3000, 'K')
            L_phot = yt.units.stefan_boltzmann_constant * 4*np.pi*radius.in_units('cm')**2 * temperature**4
            
            #Find accretion rates
            usable_inds = []
            for system in curr_components:
                if system in prev_components and system in next_components:
                    prev_ind = prev_components.index(system)
                    next_ind = next_components.index(system)
                    dM = yt.YTQuantity(next_mass[next_ind] - prev_mass[prev_ind], 'Msun')
                    dt = next_time - prev_time
                    M_dot = dM/dt
                    # This first cut is to try and pick stars that are still in protostellar stages (i.e. still accreting)
                    f_acc = 1 #Assumed all accretion energy is radiated away
                    curr_ind = curr_components.index(system)
                    M = yt.YTQuantity(curr_mass[curr_ind], 'Msun')
                    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
                    L_tot = (L_phot + L_acc.in_units('erg/s'))/yt.units.lsun.in_units('erg/s')
                    if M_dot > 1.e-7 and L_tot > 0.04 and L_tot < 32.0:
                        usable_inds.append(curr_inds[curr_components.index(system)])
                    #Calculate accretion luminosity
                    
            
            #Now bin data and create plots:
            bins = np.logspace(1.25,4,12)
            if args.projected_separation == 'False':
                current_separations = curr_res['separation'][usable_inds]
                xlabel = 'Separation (Log(AU))'
            else:
                xlabel = 'Projected Separation (Log(AU))'
                if args.axis == 'x':
                    axis_ind = 0
                elif args.axis == 'y':
                    axis_ind = 1
                else:
                    axis_ind = 2
                current_separations = np.sqrt(curr_res['separation']**2 - ((curr_res['abspos'][curr_res['index1']] - curr_res['abspos'][curr_res['index2']])[:,axis_ind])**2)
                current_separations = current_separations[usable_inds]
            cf_array = []
            n_systems = []
            s = np.where(curr_res['n']==1)[0]
            ns = np.count_nonzero(curr_res['topSystem'][s])
            for bin_it in range(1,len(bins)):
                bin_inds = np.where((current_separations>=bins[bin_it-1])&(current_separations<bins[bin_it]))[0]
                
                b = np.where(curr_res['n'][usable_inds][bin_inds]==2)[0]
                nb = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][b])
                
                t = np.where(curr_res['n'][usable_inds][bin_inds]==3)[0]
                nt = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][t])
                
                q = np.where(curr_res['n'][usable_inds][bin_inds]==4)[0]
                nq = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][q])
                
                q5 = np.where(curr_res['n'][usable_inds][bin_inds]==5)[0]
                nq5 = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][q5])
                
                s6 = np.where(curr_res['n'][usable_inds][bin_inds]==6)[0]
                ns6 = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][s6])
                
                h = np.where(curr_res['n'][usable_inds][bin_inds]==7)[0]
                nh = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][h])
                
                o = np.where(curr_res['n'][usable_inds][bin_inds]==8)[0]
                no = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][o])
                
                higher_order = np.where(curr_res['n'][usable_inds][bin_inds]>8)[0]
                nhigher = np.count_nonzero(curr_res['topSystem'][usable_inds][bin_inds][higher_order])
                
                cf = (nb+nt*2+nq*3+nq5*4+ns6*5+nh*6+no*7+nhigher*8)/(ns+nb+nt+nq+nq5+ns6+nh+no+nhigher)
                n_systems.append([ns,nb,nt,nq,nq5,ns6,nh,no,nhigher])
                cf_array.append(cf)
                
            Separations.append(current_separations)
            Times.append(curr_time)
            
            #Create histograms:
            plt.clf()
            plt.hist(Separations, bins=np.log10(bins))
            plt.xlabel(xlabel)
            plt.ylabel('Frequency')
            plt.xlim([1,4])
            plt.savefig(savedir + 'separation_histogram_' + str(int(np.round(curr_time))) + '.jpg')
            print('created separation_histogram_' + str(int(np.round(curr_time))) + '.jpg')
            
            plt.clf()
            plt.bar(((np.log10(bins[:-1])+np.log10(bins[1:]))/2), cf_array, width=0.25, fill=False, edgecolor='black')
            plt.xlabel(xlabel)
            plt.ylabel('Companion Frequency')
            plt.xlim([1,4])
            plt.ylim([0.0, 0.25])
            plt.text(1.1, 0.23, 't='+str(int(np.round(curr_time)))+'$yr$')
            plt.savefig(savedir + 'companion_frequency_' + str(int(np.round(curr_time))) + '.jpg')
            print('created companion_frequency_' + str(int(np.round(curr_time))) + '.jpg')
            CF_Array_Full.append(cf_array)
            N_sys_total.append(n_systems)
            
            #Update saved arrays:
            prev_time = curr_time
            prev_inds = curr_inds
            prev_components = curr_components
            prev_mass = curr_mass
            curr_time = next_time
            curr_inds = next_inds
            curr_components = next_components
            curr_mass = next_mass
            curr_res = res
    except:
        print('No sinks found')
            
            
#Print Calculate summed values
N_sys_total = np.array(N_sys_total)
N_total = np.sum(N_sys_total,axis=0)

CF_bottom = np.sum(N_total, axis=1)
CF_top = N_total[:, 1] + N_total[:, 2]*2 + N_total[:, 3]*3 + N_total[:, 4]*4 + N_total[:, 5]*5 + N_total[:, 6]*6 + N_total[:, 7]*7 + N_total[:, 8]*8
CF_total = CF_top/CF_bottom

plt.clf()
plt.bar(((np.log10(bins[:-1])+np.log10(bins[1:]))/2), CF_total, width=0.25, fill=False, edgecolor='black')
plt.xlabel(xlabel)
plt.ylabel('Companion Frequency')
plt.xlim([1,4])
plt.ylim([0.0, 0.25])
plt.savefig(savedir + 'Total_companion_frequency_.jpg')
print('created Total_companion_frequency.jpg')

CF_Array_Full = np.array(CF_Array_Full)
CF_sum = np.sum(CF_Array_Full, axis=0)

plt.clf()
plt.bar(((np.log10(bins[:-1])+np.log10(bins[1:]))/2), CF_sum, width=0.25, fill=False, edgecolor='black')
plt.xlabel(xlabel)
plt.ylabel('Companion Frequency')
plt.xlim([1,4])
plt.savefig(savedir + 'Sum_companion_frequency_.jpg')
print('created Sum_companion_frequency.jpg')

"""
else:
    for nout in range(len(files)):
        systems={'id':[],
        'components':[],
        'dx':[],
        'dy':[],
        'dz':[],
        'separation':[],
        'reduced_mass':[],
        'x':[],
        'y':[],
        'z':[],
        'ux':[],
        'uy':[],
        'uz':[]}
        
        try:
            sink_data = rsink(nout, datadir=datadir)
        except:
            print('No sink data from nout', nout)
        
        if 'sink_data' in locals():
            sink_id = 0
            while sink_id < len(sink_data['m']):
                dx = (sink_data['x'][sink_id] - sink_data['x'])*units['length_unit'].in_units('cm')
                dy = (sink_data['y'][sink_id] - sink_data['y'])*units['length_unit'].in_units('cm')
                dz = (sink_data['z'][sink_id] - sink_data['z'])*units['length_unit'].in_units('cm')
                separation = np.sqrt(dx**2 + dy**2 + dz**2)
                
                vx_com = ((sink_data['ux'][sink_id]*sink_data['m'][sink_id] + sink_data['ux']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))*units['velocity_unit'].in_units('cm/s')
                vy_com = ((sink_data['uy'][sink_id]*sink_data['m'][sink_id] + sink_data['uy']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))*units['velocity_unit'].in_units('cm/s')
                vz_com = ((sink_data['uz'][sink_id]*sink_data['m'][sink_id] + sink_data['uz']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))*units['velocity_unit'].in_units('cm/s')
                
                dvx_1 = sink_data['ux'][sink_id]*units['velocity_unit'].in_units('cm/s') - vx_com
                dvy_1 = sink_data['uy'][sink_id]*units['velocity_unit'].in_units('cm/s') - vy_com
                dvz_1 = sink_data['uz'][sink_id]*units['velocity_unit'].in_units('cm/s') - vz_com
                
                dvx_2 = sink_data['ux']*units['velocity_unit'].in_units('cm/s') - vx_com
                dvy_2 = sink_data['uy']*units['velocity_unit'].in_units('cm/s') - vy_com
                dvz_2 = sink_data['uz']*units['velocity_unit'].in_units('cm/s') - vz_com
                
                E_pot = (-1*(yt.units.G*(sink_data['m'][sink_id]*sink_data['m'])*units['mass_unit'].in_units('g')**2)/separation.in_units('cm')).in_units('erg')
                E_kin = (0.5*sink_data['m'][sink_id]*units['mass_unit'].in_units('g')*(dvx_1**2 + dvy_1**2 + dvz_1**2) + 0.5*sink_data['m']*units['mass_unit'].in_units('g')*(dvx_2**2 + dvy_2**2 + dvz_2**2)).in_units('erg')
        
                E_tot = E_pot + E_kin
                bound_inds = np.where((E_tot<0) & (np.isinf(E_tot)==False))[0]
                if sink_id in systems['id']:
                    sys_it = systems['id'].index(sink_id)
                    for comp in systems['components'][sys_it]:
                        if comp in bound_inds:
                            del_ind = np.argwhere(bound_inds == comp)
                            bound_inds = np.delete(bound_inds, del_ind)
                if 22 in bound_inds:
                    import pdb
                    pdb.set_trace()

                if len(bound_inds) > 0:
                    most_bound_ind = np.where(E_tot == np.min(E_tot[bound_inds]))[0][0]
                    dx = (sink_data['x'][most_bound_ind] - sink_data['x'])*units['length_unit'].in_units('cm')
                    dy = (sink_data['y'][most_bound_ind] - sink_data['y'])*units['length_unit'].in_units('cm')
                    dz = (sink_data['z'][most_bound_ind] - sink_data['z'])*units['length_unit'].in_units('cm')
                    separation = np.sqrt(dx**2 + dy**2 + dz**2)
                    
                    vx_com = ((sink_data['ux'][most_bound_ind]*sink_data['m'][most_bound_ind] + sink_data['ux']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))
                    vy_com = ((sink_data['uy'][most_bound_ind]*sink_data['m'][most_bound_ind] + sink_data['uy']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))
                    vz_com = ((sink_data['uz'][most_bound_ind]*sink_data['m'][most_bound_ind] + sink_data['uz']*sink_data['m'])/(sink_data['m'][sink_id]+sink_data['m']))
                    
                    dvx_1 = (sink_data['ux'][most_bound_ind] - vx_com)*units['velocity_unit'].in_units('cm/s')
                    dvy_1 = (sink_data['uy'][most_bound_ind] - vy_com)*units['velocity_unit'].in_units('cm/s')
                    dvz_1 = (sink_data['uz'][most_bound_ind] - vz_com)*units['velocity_unit'].in_units('cm/s')
                    
                    dvx_2 = (sink_data['ux'] - vx_com)*units['velocity_unit'].in_units('cm/s')
                    dvy_2 = (sink_data['uy'] - vy_com)*units['velocity_unit'].in_units('cm/s')
                    dvz_2 = (sink_data['uz'] - vz_com)*units['velocity_unit'].in_units('cm/s')
                    
                    E_pot = (-1*(yt.units.G*(sink_data['m'][most_bound_ind]*sink_data['m'])*units['mass_unit'].in_units('g')**2)/separation.in_units('cm')).in_units('erg')
                    E_kin = (0.5*sink_data['m'][most_bound_ind]*units['mass_unit'].in_units('g')*(dvx_1**2 + dvy_1**2 + dvz_1**2) + 0.5*sink_data['m']*units['mass_unit'].in_units('g')*(dvx_2**2 + dvy_2**2 + dvz_2**2)).in_units('erg')
            
                    E_tot = E_pot + E_kin
                    if np.argsort(E_tot)[1] == sink_id and E_tot[sink_id]<0:
                        print('Found Binary')
                        
                        system_id = len(sink_data['m'])
                        systems['id'].append(system_id)
                        systems['components'].append([sink_id, most_bound_ind])
                        systems['dx'].append(dx[sink_id])
                        systems['dy'].append(dy[sink_id])
                        systems['dz'].append(dz[sink_id])
                        systems['separation'].append(separation[sink_id].in_units('AU'))
                        
                        reduced_mass = np.product(sink_data['m'][[sink_id,most_bound_ind]])/np.sum(sink_data['m'][[sink_id,most_bound_ind]])
                        systems['reduced_mass'].append(reduced_mass)
                        sink_data['m'] = np.concatenate((sink_data['m'],[reduced_mass]))
                        x_com = (np.sum(sink_data['x'][[sink_id,most_bound_ind]]*sink_data['m'][[sink_id,most_bound_ind]])/np.sum(sink_data['m'][[sink_id,most_bound_ind]]))
                        y_com = (np.sum(sink_data['y'][[sink_id,most_bound_ind]]*sink_data['m'][[sink_id,most_bound_ind]])/np.sum(sink_data['m'][[sink_id,most_bound_ind]]))
                        z_com = (np.sum(sink_data['z'][[sink_id,most_bound_ind]]*sink_data['m'][[sink_id,most_bound_ind]])/np.sum(sink_data['m'][[sink_id,most_bound_ind]]))
                        
                        systems['x'].append(x_com)
                        sink_data['x'] = np.concatenate((sink_data['x'],[x_com]))
                        systems['y'].append(y_com)
                        sink_data['y'] = np.concatenate((sink_data['y'],[y_com]))
                        systems['z'].append(z_com)
                        sink_data['z'] = np.concatenate((sink_data['z'],[z_com]))
                        systems['ux'].append(vx_com[sink_id])
                        sink_data['ux'] = np.concatenate((sink_data['ux'],[vx_com[sink_id]]))
                        systems['uy'].append(vy_com[sink_id])
                        sink_data['uy'] = np.concatenate((sink_data['uy'],[vy_com[sink_id]]))
                        systems['uz'].append(vz_com[sink_id])
                        sink_data['uz'] = np.concatenate((sink_data['uz'],[vz_com[sink_id]]))
                else:
                    print("no bound sinks found for sink", sink_id)
                sink_id = sink_id + 1
                '''
                if len(systems['id'])>1:
                    import pdb
                    pdb.set_trace()
                '''
"""
