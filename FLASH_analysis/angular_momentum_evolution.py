import numpy as np
import pickle
import yt
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
from mpi4py.MPI import COMM_WORLD as CW
import my_flash_fields as myf
import my_flash_module as mym

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

line_styles = ['--', '-.', '-']
label = ['Primary', 'Secondary', 'Orbit']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

#---------------------------------------------------

#Get simulation files
input_dir = sys.argv[1]
files = sorted(glob.glob(input_dir + '*plt_cnt*'))

#Get first file with sink:
start_file = mym.find_files([0], files)[0]
files = files[files.index(start_file):files.index(start_file)+10] #files[files.index(start_file):]
ts = yt.DatasetSeries(files, parallel=True)

L_dict = {}
Time_array = []
L_primary = []
L_secondary = []
L_orbit = []
L_in_gas = []
for sto, ds in ts.piter(storage=L_dict):
    Time_array.append(ds.current_time.in_units('yr'))

    #Calculate CoM
    dd = ds.all_data()
    
    #Calculate particle spin
    particle_spin = dd['Particle_Spin']
    
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
    
    rank_data = {'Time_array': Time_array, 'L_primary': L_primary, 'L_secondary': L_secondary, 'L_orbit:': L_orbit, 'L_in_gas': L_in_gas}
    sto.result_id = 'rank_'+str(rank)
    sto.result = rank_data
    

#Compile together results
Time_array = []
L_primary = []
L_secondary = []
L_orbit = []
L_in_gas = []
for key in L_dict.keys():
    Time_array = Time_array + L_dict[key]['Time_array']
    L_primary = L_primary + L_dict[key]['L_primary']
    L_secondary = L_secondary + L_dict[key]['L_secondary']
    L_orbit = L_orbit + L_dict[key]['L_orbit']
    L_in_gas = L_in_gas + L_dict[key]['L_in_gas']

#sort arrays
sorted_inds = np.argsort(Time_array)
Time_array = np.array(Time_array)[sorted_inds]
L_primary = np.array(L_primary)[sorted_inds]
L_secondary = np.array(L_secondary)[sorted_inds]
L_orbit = np.array(L_orbit)[sorted_inds]
L_in_gas = np.array(L_in_gas)[sorted_inds]

plt.clf()
plt.semilogy(Time_array, L_primary, label='Primary spin')
plt.semilogy(Time_array, L_secondary, label='Secondary spin')
plt.semilogy(Time_array, L_orbit, label='Orbital L')
plt.semilogy(Time_array, L_in_gas, label='L_in_gas')
plt.xlabel('Time (yr)')
plt.ylabel('Angular momentum (g cm$^2$/s)')
plt.legend(loc='best')
plt.savefig('L_evolution.png', bbox_inches='tight')
