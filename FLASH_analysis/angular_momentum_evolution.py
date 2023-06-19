import numpy as np
import pickle
import yt
yt.enable_parallelism()
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
from mpi4py.MPI import COMM_WORLD as CW
import my_flash_fields as myf
import my_flash_module as mym
import pickle
import argparse
import os

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
#Ploting parameters
matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

#---------------------------------------------------
#Define arguments
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-update", "--update_pickles", help="do you want to read the Flash output and update the pickles?", type=str, default='True')
    parser.add_argument("-dt", "--time_step", help="what time step between data points do you want to use?", type=float, default=10.0)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
input_dir = sys.argv[1]
args = parse_inputs()
files = sorted(glob.glob(input_dir + '*plt_cnt*'))

if args.update_pickles == 'True':

    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None)

    L_dict = {}
    Time_array = []
    L_sink = {}
    L_orbit = []
    L_in_gas = []
    T_round_all = []

    pickle_names = '_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'ang_mom_*.pkl'
    #get current progress
    if rank == 0:
        pickle_files = glob.glob(pickle_names)
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                L_dict = pickle.load(file)
                file.close()
                
                for key in L_dict.keys():
                    Time_array = Time_array + L_dict['Time_array']
                    import pdb
                    pdb.set_trace()
                    L_orbit = L_orbit + L_dict['L_orbit']
                    L_in_gas = L_in_gas + L_dict['L_in_gas']
                    T_round_all = T_round_all + L_dict['T_round']
            
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            import pdb
            pdb.set_trace()
            L_orbit = list(np.array(L_orbit)[sorted_inds])
            L_in_gas = list(np.array(L_in_gas)[sorted_inds])
            
            start_time = np.max(T_round_all)
        else:
            start_time = m_times[0]
    else:
        start_time = np.nan

    sys.stdout.flush()
    CW.Barrier()

    start_time = CW.bcast(start_time, root=0)

    sys.stdout.flush()
    CW.Barrier()
    
    no_frames = len(m_times)
    start_frame = m_times.index(start_time)
    m_times = m_times[start_frame:]
    usable_files = mym.find_files(m_times, files)
    frames = list(range(start_frame, no_frames))

    #make time series
    ts = yt.DatasetSeries(usable_files, parallel=True)
    form_time = mym.find_sink_formation_time(files)
    form_time = yt.YTQuantity(form_time, 'yr')

    sys.stdout.flush()
    CW.Barrier()

    for ds in ts.piter():
        t_round = m_times[usable_files.index(ds.filename)]
        time_val = ds.current_time.in_units('yr') - form_time
        Time_array.append(time_val)

        #load all data
        dd = ds.all_data()
        
        if time_val > 0:
            #Calculate particle spin
            #Define sink pickle
            sink_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_'+pickle_names.split('_ang_mom_')[0] + '.pkl'
            file = open(sink_pickle, 'rb')
            sink_data = pickle.load(file)
            file.close()
            
            prime_tag = list(sink_data.keys())[0]
            match_time_ind = np.argmin(abs(sink_data[prime_tag]['time']- ds.current_time.value))
            prime_spin = np.sqrt(sink_data[prime_tag]['anglx'][match_time_ind]**2 + sink_data[prime_tag]['angly'][match_time_ind]**2 + sink_data[prime_tag]['anglz'][match_time_ind]**2)
            
            try:
                particle_spin = yt.YTArray(np.sqrt(dd['particle_x_ang']**2 + dd['particle_y_ang']**2 + dd['particle_z_ang']**2), 'g*cm**2/s')
            except:
                particle_spin = yt.YTArray([np.sqrt(dd['particle_x_ang']**2 + dd['particle_y_ang']**2 + dd['particle_z_ang']**2)], 'g*cm**2/s')
            
            particle_tags = dd['particle_tag']
            tag_sort_inds = []
            for sink_tag in list(sink_data.keys()):
                if int(sink_tag) in particle_tags:
                    tag_sort_inds.append(list(particle_tags.value).index(float(sink_tag)))
            
            particle_spin = particle_spin[tag_sort_inds]
            
            if '00000' in str(prime_spin):
                spin_string = ''.join(str(prime_spin).split('00000'))
            else:
                spin_string = str(prime_spin)
            round_int = len(spin_string.split('.')[-1].split('e')[0])
            format_string_line = str("'{:0."+str(round_int)+"e}'.format(particle_spin[0].value)")
            if eval(format_string_line) != spin_string:
                if abs(float(eval(format_string_line)) - float(spin_string)) > 1.e45:
                    dt_discrep = yt.YTQuantity(abs(ds.current_time.value - sink_data[prime_tag]['time'][match_time_ind]), 's').in_units('yr')
                    if dt_discrep < 5:
                        print("Spin is diverging between YT and the Sink_evol.dat!")
                        #if size == 1:
                        #    import pdb
                        #    pdb.set_trace()
                
            
            #Calculate orbital angular momentum around CoM
            dx = dd['particle_posx'].in_units('cm') - dd['CoM_full'][0]
            dy = dd['particle_posy'].in_units('cm') - dd['CoM_full'][1]
            dz = dd['particle_posz'].in_units('cm') - dd['CoM_full'][2]
            d_pos = yt.YTArray([dx, dy, dz]).T
            
            dvx = dd['particle_velx'].in_units('cm/s') - dd['CoM_Velocity_full'][0]
            dvy = dd['particle_vely'].in_units('cm/s') - dd['CoM_Velocity_full'][1]
            dvz = dd['particle_velz'].in_units('cm/s') - dd['CoM_Velocity_full'][2]
            d_vel = yt.YTArray([dvx, dvy, dvz]).T
            
            L_orb = dd['particle_mass'].value * np.cross(d_vel, d_pos).T
            L_orb_tot = yt.YTQuantity(np.sum(np.sqrt(np.sum(L_orb**2, axis=0))), 'g*cm**2/s')
        else:
            particle_spin = yt.YTArray([np.nan], 'g*cm**2/s')
            L_orb_tot = yt.YTArray(np.nan, 'g*cm**2/s')
        
        #Calculate angular momentum in gas
        dx_gas = dd['x'] - dd['CoM_full'][0]
        dy_gas = dd['y'] - dd['CoM_full'][1]
        dz_gas = dd['z'] - dd['CoM_full'][2]
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        
        dvx_gas = dd['velx'].in_units('cm/s') - dd['CoM_Velocity_full'][0]
        dvy_gas = dd['vely'].in_units('cm/s') - dd['CoM_Velocity_full'][1]
        dvz_gas = dd['velz'].in_units('cm/s') - dd['CoM_Velocity_full'][2]
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        
        L_gas = dd['mass'].value * np.cross(d_vel_gas, d_pos_gas).T
        L_gas_tot = yt.YTQuantity(np.sum(np.sqrt(np.sum(L_gas**2, axis=0))), 'g*cm**2/s')
        
        #Save values
        for particle_tag in particle_tags:
            if particle_tag not in L_sink.keys():
                L_sink.update({str(particle_tag):[time_val, particle_spin[list(particle_tags).index(particle_tag)]]})
            else:
                L_sink[str(particle_tag)].append([time_val, particle_spin[list(particle_tags).index(particle_tag)]])
        '''
        L_primary.append(particle_spin[0])
        if len(particle_spin) == 2:
            L_secondary.append(particle_spin[1])
            
            sec_tag = list(sink_data.keys())[1]
            match_time_ind = np.argmin(abs(sink_data[sec_tag]['time']- ds.current_time.value))
            sec_spin = np.sqrt(sink_data[sec_tag]['anglx'][match_time_ind]**2 + sink_data[sec_tag]['angly'][match_time_ind]**2 + sink_data[sec_tag]['anglz'][match_time_ind]**2)
            
            if '00000' in str(sec_spin):
                spin_string = ''.join(str(sec_spin).split('00000'))
            else:
                spin_string = str(sec_spin)
            round_int = len(spin_string.split('.')[-1].split('e')[0])
            format_string_line = str("'{:0."+str(round_int)+"e}'.format(particle_spin[1].value)")
            if eval(format_string_line) != spin_string:
                if abs(float(eval(format_string_line)) - float(spin_string)) > 1.e45:
                    dt_discrep = yt.YTQuantity(abs(ds.current_time.value - sink_data[prime_tag]['time'][match_time_ind]), 's').in_units('yr')
                    if dt_discrep < 5:
                        print("Spin is diverging between YT and the Sink_evol.dat!")
                        #if size == 1:
                        #    import pdb
                        #    pdb.set_trace()
        else:
            L_secondary.append(yt.YTQuantity(np.nan, particle_spin.units))
        '''
        L_orbit.append(L_orb_tot)
        L_in_gas.append(L_gas_tot)
        
        rank_data = {'Time_array': Time_array, 'L_orbit': L_orbit, 'L_in_gas': L_in_gas, 'T_round': [t_round], 'L_sink': L_sink}
        
        #write pickle
        file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'ang_mom_'+str(rank)+'.pkl', 'wb')
        pickle.dump((rank_data), file)
        file.close()
        print('read file', usable_files.index(ds.filename), 'of', len(usable_files), 'files on rank', rank)

    sys.stdout.flush()
    CW.Barrier()

    if rank == 0:
        #Compile together results
        Time_array_full = []
        L_orbit_full = []
        L_in_gas_full = []
        L_sink_full = {}
        pickle_files = glob.glob('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'ang_mom_*.pkl')
        for pickle_file in pickle_files:
            file = open(pickle_file, 'rb')
            rank_data = pickle.load(file)
            file.close()
            
            for key in rank_data.keys():
                Time_array_full = Time_array_full + rank_data['Time_array']
                L_orbit_full = L_orbit_full + rank_data['L_orbit']
                L_in_gas_full = L_in_gas_full + rank_data['L_in_gas']
                for key in rank_data['L_sink'].keys():
                    if key not in L_sink_full.keys():
                        L_sink_full.update({key:rank_data['L_sink'][key]})
                    else:
                        L_sink_full[key].append(rank_data['L_sink'][key])
        
        sorted_inds = np.argsort(Time_array_full)
        Time_array = np.array(Time_array_full)[sorted_inds]
        for key in L_sink_full:
            t_sorted_inds = np.argsort(L_sink_full[key].T[0])
            L_sink_full[key] = L_sink_full[key][t_sorted_inds]
        L_orbit = np.array(L_orbit)[sorted_inds]
        L_in_gas = np.array(L_in_gas)[sorted_inds]
        
        file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'gathered_ang_mom.pkl', 'wb')
        pickle.dump((Time_array, L_sink_full, L_orbit, L_in_gas), file)
        file.close()
        
        for pickle_file in pickle_files:
            os.remove(pickle_file)
        print('saved gathered data')

sys.stdout.flush()
CW.Barrier()

if rank == 0:
    file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'gathered_ang_mom.pkl', 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.clf()
    #plt.semilogy(Time_array - Time_array[0], L_primary, label='Primary spin')
    #plt.semilogy(Time_array - Time_array[0], L_secondary, label='Secondary spin')
    plt.semilogy(Time_array - Time_array[0], L_orbit, label='Orbital L')
    plt.semilogy(Time_array - Time_array[0], L_in_gas, label='L_in_gas')
    plt.xlabel('Time (yr)')
    plt.ylabel('Angular momentum (g cm$^2$/s)')
    plt.legend(loc='best')
    plt.xlim(left=0)
    plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
    plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.savefig('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'L_evolution.png', bbox_inches='tight')
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    plt.clf()
    #plt.semilogy(Time_array - Time_array[0], L_primary/L_tot, label='Primary spin')
    #plt.semilogy(Time_array - Time_array[0], L_secondary/L_tot, label='Secondary spin')
    plt.semilogy(Time_array - Time_array[0], L_orbit/L_tot, label='Orbital L')
    plt.semilogy(Time_array - Time_array[0], L_in_gas/L_tot, label='L_in_gas')
    plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
    plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.xlabel('Time (yr)')
    plt.ylabel('Angular momentum fraction (%)')
    plt.legend(loc='best')
    plt.xlim(left=0)
    plt.ylim(top=1)
    plt.savefig('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'L_evolution_frac.png', bbox_inches='tight')
    
    
    print('saved figure')

sys.stdout.flush()
CW.Barrier()
