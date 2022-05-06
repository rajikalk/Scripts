import numpy as np
import matplotlib.pyplot as plt
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

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

#Formation_pathway = [[252, 88, 105], [624, 859, 1209], [1893, 2000, 5105], [1255, 4381, 11458], [2921, 5951, 22602], [2172, 5412, 32412]]
birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/sink_birth_all.pkl"]


Core_frag_fracs = []
Delayed_core_frag_fracs = []
Dynamical_capt_fracs = []
Initial_Seps = []
for birth_con_pickle in birth_con_pickles:
    file = open(birth_con_pickle, 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    total_components = 0
    Core_frag_counter = 0
    Delayed_core_fra_counter = 0
    Dynamical_capt_counter = 0
    for key in Sink_birth_all.keys():
        total_components = total_components + 1
        if Sink_birth_all[str(np.max(sub_sys))][0] == True:
            Core_frag_counter = Core_frag_counter + 1
        elif str(Sink_birth_all[str(np.max(sub_sys))][1]) in Sink_birth_all[str(np.max(sub_sys))][2]:
            Delayed_core_fra_counter = Delayed_core_fra_counter + 1
        else:
            Dynamical_capt_counter = Dynamical_capt_counter + 1
        
    Core_frag_frac = Core_frag_counter/total_components
    Delayed_core_frag_frac = Delayed_core_frag_frac/total_components
    Dynamical_capt_frac = Dynamical_capt_frac/total_components
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)

x_labels = ['Core frag.', 'Delayed core frag.', 'Dynamical capture']

fig, ax = plt.subplots(1, 1, figsize=(5, 4))

ind = np.arange(len(subplot_titles))
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

p1 = plt.bar(ind, Core_frag_fracs, 0.95, color='b', linewidth=1, edgecolor='k')#, hatch='+'
p2 = plt.bar(ind, Delayed_core_frag_fracs, 0.95, bottom=Core_frag_fracs, color='m', linewidth=1, edgecolor='k')#, hatch='x'
p3 = plt.bar(ind, Dynamical_capt_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Core_frag_fracs)), color='r', linewidth=1, edgecolor='k')#, hatch='O'

plt.xlim([-0.6, 5.6])
plt.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=font_size, right=True)
ax.tick_params(axis='both', which='minor', labelsize=font_size, left=True, right=True, top=False, bottom=False)
plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
ax.tick_params(which='both', direction='in')
plt.xlabel('Initial Gas Mass (M$_\odot$)', fontsize=font_size, labelpad=-0.5)
plt.legend((p3[0], p2[0], p1[0]), ('Dynamical capture', 'Delayed core frag.', 'Core fragmentation'), loc='upper right', fontsize=font_size)
plt.ylabel('Fraction', fontsize=font_size, labelpad=-0.5)
plt.ylim([0,1])
plt.savefig('formation_pathway.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

