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
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
input_dir = sys.argv[1]
args = parse_inputs()
files = sorted(glob.glob(input_dir + '*plt_cnt*'))

if args.update_pickles == 'True':

    L_dict = {}
    Time_array = []
    L_primary = []
    L_secondary = []
    L_orbit = []
    L_in_gas = []

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
                    L_primary = L_primary + L_dict['L_primary']
                    L_secondary = L_secondary + L_dict['L_secondary']
                    L_orbit = L_orbit + L_dict['L_orbit']
                    L_in_gas = L_in_gas + L_dict['L_in_gas']
            
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            L_primary = list(np.array(L_primary)[sorted_inds])
            L_secondary = list(np.array(L_secondary)[sorted_inds])
            L_orbit = list(np.array(L_orbit)[sorted_inds])
            L_in_gas = list(np.array(L_in_gas)[sorted_inds])
            
            start_time = Time_array[-1] - mym.find_sink_formation_time(files)
        else:
            start_time = 0
        start_file = mym.find_files([start_time], files)[0]
    else:
        start_file = ''

    sys.stdout.flush()
    CW.Barrier()

    start_file = CW.bcast(start_file, root=0)

    sys.stdout.flush()
    CW.Barrier()

    #make time series
    files = files[files.index(start_file):] #files[files.index(start_file):]
    ts = yt.DatasetSeries(files, parallel=True)

    sys.stdout.flush()
    CW.Barrier()

    for ds in ts.piter():
        Time_array.append(ds.current_time.in_units('yr'))

        #load all data
        dd = ds.all_data()
        
        #Calculate particle spin
        import pdb
        pdb.set_trace()
        particle_spin = dd['particle_angular_momentum_magnitude']
        
        #Calculate orbital angular momentum around CoM
        dx = dd['particle_posx'].in_units('cm') - dd['CoM'][0]
        dy = dd['particle_posy'].in_units('cm') - dd['CoM'][1]
        dz = dd['particle_posz'].in_units('cm') - dd['CoM'][2]
        d_pos = yt.YTArray([dx, dy, dz]).T
        
        dvx = dd['particle_velx'].in_units('cm/s') - dd['CoM_Velocity'][0]
        dvy = dd['particle_vely'].in_units('cm/s') - dd['CoM_Velocity'][1]
        dvz = dd['particle_velz'].in_units('cm/s') - dd['CoM_Velocity'][2]
        d_vel = yt.YTArray([dvx, dvy, dvz]).T
        
        L_orb = dd['particle_mass'].value * np.cross(d_vel, d_pos).T
        L_orb_tot = yt.YTQuantity(np.sum(np.sqrt(np.sum(L_orb**2, axis=0))), 'g*cm**2/s')
        
        #Calculate angular momentum in gas
        dx_gas = dd['x'] - dd['CoM'][0]
        dy_gas = dd['y'] - dd['CoM'][1]
        dz_gas = dd['z'] - dd['CoM'][2]
        d_pos_gas = yt.YTArray([dx_gas, dy_gas, dz_gas]).T
        
        dvx_gas = dd['velx'].in_units('cm/s') - dd['CoM_Velocity'][0]
        dvy_gas = dd['vely'].in_units('cm/s') - dd['CoM_Velocity'][1]
        dvz_gas = dd['velz'].in_units('cm/s') - dd['CoM_Velocity'][2]
        d_vel_gas = yt.YTArray([dvx_gas, dvy_gas, dvz_gas]).T
        
        L_gas = dd['mass'].value * np.cross(d_vel_gas, d_pos_gas).T
        L_gas_tot = yt.YTQuantity(np.sum(np.sqrt(np.sum(L_gas**2, axis=0))), 'g*cm**2/s')
        
        #Save values
        L_primary.append(particle_spin[0])
        if len(particle_spin) == 2:
            L_secondary.append(particle_spin[1])
        else:
            L_secondary.append(yt.YTQuantity(np.nan, particle_spin.units))
        L_orbit.append(L_orb_tot)
        L_in_gas.append(L_gas_tot)
        
        rank_data = {'Time_array': Time_array, 'L_primary': L_primary, 'L_secondary': L_secondary, 'L_orbit': L_orbit, 'L_in_gas': L_in_gas}
        
        #write pickle
        file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'ang_mom_'+str(rank)+'.pkl', 'wb')
        pickle.dump((rank_data), file)
        file.close()
        print('read file', files.index(ds.filename), 'of', len(files), 'files on rank', rank)

    sys.stdout.flush()
    CW.Barrier()

    if rank == 0:
        #Compile together results
        pickle_files = glob.glob('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'ang_mom_*.pkl')
        for pickle_file in pickle_files:
            file = open(pickle_file, 'rb')
            rank_data = pickle.load(file)
            file.close()
            
            for key in rank_data.keys():
                Time_array = Time_array + rank_data['Time_array']
                L_primary = L_primary + rank_data['L_primary']
                L_secondary = L_secondary + rank_data['L_secondary']
                L_orbit = L_orbit + rank_data['L_orbit']
                L_in_gas = L_in_gas + rank_data['L_in_gas']
        
        sorted_inds = np.argsort(Time_array)
        Time_array = np.array(Time_array)[sorted_inds]
        L_primary = np.array(L_primary)[sorted_inds]
        L_secondary = np.array(L_secondary)[sorted_inds]
        L_orbit = np.array(L_orbit)[sorted_inds]
        L_in_gas = np.array(L_in_gas)[sorted_inds]
        
        file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'gathered_ang_mom.pkl', 'wb')
        pickle.dump((Time_array, L_primary, L_secondary, L_orbit, L_in_gas), file)
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
    plt.semilogy(Time_array - Time_array[0], L_primary, label='Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary, label='Secondary spin')
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
    plt.semilogy(Time_array - Time_array[0], L_primary/L_tot, label='Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary/L_tot, label='Secondary spin')
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
