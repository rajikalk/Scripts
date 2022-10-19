import numpy as np
import pickle
import yt
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
#from mpi4py.MPI import COMM_WORLD as CW

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

input_dir = sys.argv[1]
pickle_files = sorted(glob.glob(input_dir+'*/Lref_11.pkl'))

line_styles = ['--', '-.', '-']
#label = ['L_{ref}=10', 'L_{ref}=11', 'L_{ref}=12']
label = ['$\Omega\times t_{ff}$=0.25, $\alpha$=0.75, $\mathcal(M)$=0.1', '$\Omega\times t_{ff}$=0.2, $\alpha$=0.50, $\mathcal(M)$=0.0', '$\Omega\times t_{ff}$=0.2, $\alpha$=0.75, $\mathcal(M)$=0.1']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()

pick_it = -1
for pickle_file in pickle_files:
    pick_it = pick_it + 1
    file_open = open(pickle_file, 'rb')
    sink_data = pickle.load(file_open)
    file_open.close()
    
    secondary_sink_form_time = sink_data[list(sink_data.keys())[1]]['time'][0]
    primary_time_star_ind = np.argwhere(sink_data[list(sink_data.keys())[0]]['time']==secondary_sink_form_time)[0][0]
    
    dx = sink_data[list(sink_data.keys())[0]]['posx'][primary_time_star_ind:] - sink_data[list(sink_data.keys())[1]]['posx']
    dy = sink_data[list(sink_data.keys())[0]]['posy'][primary_time_star_ind:] - sink_data[list(sink_data.keys())[1]]['posy']
    dz = sink_data[list(sink_data.keys())[0]]['posz'][primary_time_star_ind:] - sink_data[list(sink_data.keys())[1]]['posz']

    separation = yt.YTArray(np.sqrt(dx**2 + dy**2 + dz**2), 'cm')
    time = sink_data[list(sink_data.keys())[1]]['time'] - sink_data[list(sink_data.keys())[1]]['time'][0]
    
    plt.plot(time.in_units('yr'), separation.in_units('au'), label=label[pick_it])

plt.xlabel('Time (yr)')
plt.ylabel('Separation (au)')
plt.legend()
plt.savefig('separation_evolution.png')
