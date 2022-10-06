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
pickle_files = sorted(glob.glob(input_dir+'*.pkl'))

line_styles = ['-.', '--', '-']
label = ['L_{ref}=10', 'L_{ref}=11', 'L_{ref}=12']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width, single_col_width*1.5), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

pick_it = -1
for pickle_file in pickle_files:
    pick_it = pick_it + 1
    file_open = open(pickle_file, 'rb')
    sink_data = pickle.load(file_open)
    file_open.close()
    
    plot_counter = 0
    first_sink_formation = np.nan
    for sink_id in sink_data.keys():
        if np.isnan(first_sink_formation):
            first_sink_formation = sink_data[sink_id]['time'][0]
        time_arr = yt.YTArray(sink_data[sink_id]['time']-first_sink_formation, 's')
        Lx = yt.YTArray(sink_data[sink_id]['anglx'], 'g*cm/s')
        Ly = yt.YTArray(sink_data[sink_id]['angly'], 'g*cm/s')
        Lz = yt.YTArray(sink_data[sink_id]['anglz'], 'g*cm/s')
        L_tot = np.sqrt(Lx**2 + Ly**2 + Lz**2)
        
        axs.flatten()[plot_counter].semilogy(time_arr.in_units('yr'), L_tot, ls=line_styles[pick_it], label=label[pick_it])
        
        axs.flatten()[plot_counter].set_ylabel('L$_{'+str(plot_counter+1)+'}$ (g cm$^2$/s)')
        plot_counter = plot_counter +1
    
    axs.flatten()[0].legend(loc='lower right')
    axs.flatten()[plot_counter-1].set_xlabel('Time (yr)')
    axs.flatten()[plot_counter-1].set_xlim(left=0)
    
plt.savefig('spin_evolution.pdf', bbox_inches='tight', pad_inches=0.02)
