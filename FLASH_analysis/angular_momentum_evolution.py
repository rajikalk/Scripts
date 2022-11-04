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
files = files[files.index(start_file):]

L_dict = {}
L_primary = []
L_secondary = []
L_orbit = []
L_gas = []
for fn in yt.parallel_objects(files, njobs=size, storage=L_dict):
    part_file = 'part'.join(fn[-1].split('plt_cnt'))
    ds = yt.load(fn[-1], particle_filename=part_file)

    #Calculate CoM
    dd = ds.all_data()
    CoM =
    
    particle_spin = dd['Particle_Spin']
    
    #Calculate orbital angular momentum around CoM
    dx = dd['particle_posx'].in_units('cm') - dd['CoM'][0]
    dy = dd['particle_posy'].in_units('cm') - dd['CoM'][1]
    dz = dd['particle_posz'].in_units('cm') - dd['CoM'][2]
    
    dvx = dd['particle_velx'].in_units('cm/s') - dd['CoM_Velocity'][0]
    dvy = dd['particle_vely'].in_units('cm/s') - dd['CoM_Velocity'][1]
    dvz = dd['particle_velz'].in_units('cm/s') - dd['CoM_Velocity'][2]
    
    import pdb
    pdb.set_trace()
    
