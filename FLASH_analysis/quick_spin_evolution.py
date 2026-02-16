import numpy as np
import pickle
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

B_2 = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_0.20_Binary_Lref_10_Restart_Mach_0.2_From_binary_formation.pkl'
B_25 = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_0.25_Binary_Lref_10_Restart_Mach_0.2_From_binary_formation.pkl'
S_2 = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_0.20_Single_Lref_10_Mach_0.2.pkl'
S_25 = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_0.25_Single_Lref_10_Mach_0.2.pkl'

file = open(S_2, 'rb')
sink_data = pickle.load(file)
file.close()

S_2_total_spin = np.sqrt(sink_data['65714']['anglx']**2 + sink_data['65714']['angly']**2 + sink_data['65714']['anglz']**2)
S_2_time = sink_data['65714']['time'] - sink_data['65714']['time'][0]

file = open(S_25, 'rb')
sink_data = pickle.load(file)
file.close()

S_25_time = sink_data['65711']['time'] - sink_data['65711']['time'][0]
S_25_total_spin = np.sqrt(sink_data['65711']['anglx']**2 + sink_data['65711']['angly']**2 + sink_data['65711']['anglz']**2)

file = open(B_2, 'rb')
sink_data = pickle.load(file)
file.close()

B_2_time_1 = sink_data['65727']['time'] - sink_data['65727']['time'][0]
B_2_time_2 = sink_data['65779']['time'] - sink_data['65779']['time'][0]
B_2_total_spin_1 = np.sqrt(sink_data['65727']['anglx']**2 + sink_data['65727']['angly']**2 + sink_data['65727']['anglz']**2)
B_2_total_spin_2 = np.sqrt(sink_data['65779']['anglx']**2 + sink_data['65779']['angly']**2 + sink_data['65779']['anglz']**2)

file = open(B_25, 'rb')
sink_data = pickle.load(file)
file.close()

B_25_time_1 = sink_data['65730']['time'] - sink_data['65730']['time'][0]
B_25_time_2 = sink_data['65776']['time'] - sink_data['65776']['time'][0]
B_25_total_spin_1 = np.sqrt(sink_data['65730']['anglx']**2 + sink_data['65730']['angly']**2 + sink_data['65730']['anglz']**2)
B_25_total_spin_2 = np.sqrt(sink_data['65776']['anglx']**2 + sink_data['65776']['angly']**2 + sink_data['65776']['anglz']**2)

plt.clf()
plt.figure(figsize=(5, 3))
plt.semilogy(S_2_time/31557600.0, S_2_total_spin, 'b-', label='Single Low L$_{init}$')
plt.semilogy(S_25_time/31557600.0, S_25_total_spin, 'r-', label='Single High L$_{init}$')
plt.semilogy(B_2_time_1/31557600.0, B_2_total_spin_1, 'b--', label='Primary Low L$_{init}$', alpha=0.5)
plt.semilogy(B_2_time_2/31557600.0, B_2_total_spin_2, 'b-.', label='Secondary Low L$_{init}$', alpha=0.5)
plt.semilogy(B_25_time_1/31557600.0, B_25_total_spin_1, 'r--', label='Primary High L$_{init}$', alpha=0.5)
plt.semilogy(B_25_time_2/31557600.0, B_25_total_spin_2, 'r-.', label='Secondary High L$_{init}$', alpha=0.5)
plt.xlabel('Time since formation (yr)')
plt.ylabel('Sink particle spin (g cm$^2$/s)')
plt.xlim(left=0)
plt.legend(ncol=2)
plt.savefig('protostellar_spin.pdf', bbox_inches='tight', pad_inches=0.02)
