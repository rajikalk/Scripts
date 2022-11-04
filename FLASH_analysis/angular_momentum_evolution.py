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

pickle_file = sys.argv[1]

line_styles = ['--', '-.', '-']
label = ['Primary', 'Secondary', 'Orbit']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()

file_open = open(pickle_file, 'rb')
sink_data = pickle.load(file_open)
file_open.close()

plot_counter = 0
first_sink_formation = np.nan
for sink_id in sink_data.keys():
    if np.isnan(first_sink_formation):
        first_sink_formation = sink_data[sink_id]['time'][0]
    time_arr = yt.YTArray(sink_data[sink_id]['time']-first_sink_formation, 's')
    Lx = yt.YTArray(sink_data[sink_id]['anglx'], 'g*cm**2/s')
    Ly = yt.YTArray(sink_data[sink_id]['angly'], 'g*cm**2/s')
    Lz = yt.YTArray(sink_data[sink_id]['anglz'], 'g*cm**2/s')
    L_tot = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    L_tot[np.where(L_tot==0)[0]]=np.nan
    
    axs.flatten()[plot_counter].semilogy(time_arr.in_units('yr'), L_tot, ls=line_styles[plot_counter], label=label[plot_counter])
    
    #axs.flatten()[plot_counter].set_ylabel('L$_{'+str(plot_counter+1)+'}$ (cm$^2$/s)')
    plot_counter = plot_counter +1

#Calculate orbital angular momentum
import pdb
pdb.set_trace()

axs.flatten()[plot_counter-1].set_xlabel('Time since Primary formation (yr)')
axs.flatten()[plot_counter-1].set_xlim(left=0)
axs.flatten()[plot_counter-1].set_ylim(bottom=1e40)
    
plt.savefig('spin_evolution_with_single.pdf', bbox_inches='tight', pad_inches=0.02)
