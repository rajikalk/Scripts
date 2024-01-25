import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-smooth_acc", "--smooth_accretion", help="do you want to smooth the accrtion", type=str, default='True')
    parser.add_argument("-window", "--smoothing_window", help="how big do you want the smoothing window to be, in term of phase", type=float, default=0.1)
    parser.add_argument("-savename", "--save_image_name", help="What would you like the plot to be saved as?", type=str, default="binary_evolution")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

end_peri = 13
args = parse_inputs()
read_file = sys.argv[1]
labels = ['B1*', 'B1', 'B2', 'B3']
panel_tag = ['a)', 'b)', 'c)', 'd)']
linestyles = [':', '-', '-', '-']
colors = ['b', 'b', 'orange', 'g']

pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_91/High_resolution/Remade_pickles/reduced_system_data.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_91/Remade_pickles/reduced_system_data.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_49/Remade_pickles/reduced_system_data.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_164/Remade_pickles/reduced_system_data.pkl']

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(single_col_width, two_col_width), sharex=True, sharey=False)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pickle_file in pickle_files:
    file_open = open(pickle_file, 'rb')
    reduced_systems_data = pickle.load(file_open)
    file_open.close()
    
    import pdb
    pdb.set_trace()
