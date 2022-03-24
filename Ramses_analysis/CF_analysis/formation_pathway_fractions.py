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
Formation_pathway = [[96, 129, 108], [624, 859, 1209], [1536, 994, 2917], [1340, 1746, 5248], [3353, 2539, 9586], [6559, 4865, 22248]]
Formation_pathway = [[18, 14, 14, 77], [28, 35, 43, 647], [55, 33, 111, 1391], [54, 56, 135, 1974], [104, 71, 242, 3775], [308, 105, 549, 8141]]

x_labels = ['Core frag.', 'Delayed core frag.', 'Dynamical capture']

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
Core_frag_fracs = []
Delayed_core_frag_fracs = []
Dynamical_capt_fracs = []
Other_fracs = []
ind = np.arange(len(subplot_titles))
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

for sim_it in range(len(subplot_titles)):
    total_components = np.sum(Formation_pathway[sim_it])
    Core_frag_frac = Formation_pathway[sim_it][0]/total_components
    Delayed_core_frag_frac = Formation_pathway[sim_it][1]/total_components
    Dynamical_capt_frac = Formation_pathway[sim_it][2]/total_components
    Other_frac = Formation_pathway[sim_it][3]/total_components
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)
    Other_fracs.append(Other_frac)

p1 = plt.bar(ind, Other_fracs, 0.95, color='grey')
p2 = plt.bar(ind, Dynamical_capt_fracs, 0.95, bottom=Other_fracs, color='r')
p3 = plt.bar(ind, Delayed_core_frag_fracs, 0.95, bottom=(np.array(Dynamical_capt_fracs)+np.array(Other_fracs)), color='m')
p4 = plt.bar(ind, Core_frag_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Dynamical_capt_fracs)+np.array(Other_fracs)), color='b')

plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
#plt.yticks(np.arange(0, 1.1,0.1))
ax.tick_params(which='both', direction='in')
plt.xlabel('Initial Gas Mass (M$_\odot$)', labelpad=-20)
#ax.set_xticks([1,2,3])
#ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor")
#ax.set_xticklabels(x_labels)#, rotation=45, ha="right", rotation_mode="anchor")
plt.legend((p4[0], p3[0], p2[0], p1[0]), ('Core fragmentation', 'Delayed core fragmentation', 'Dynamical capture', 'Other'), loc='lower right')
plt.ylabel('Fraction', labelpad=-20)
plt.ylim([0,1])
plt.savefig('formation_pathway_other.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

fig, ax = plt.subplots(1, 1, figsize=(single_col_width, single_col_width-0.5))
Core_frag_fracs = []
Delayed_core_frag_fracs = []
Dynamical_capt_fracs = []
ind = np.arange(len(subplot_titles))

for sim_it in range(len(subplot_titles)):
    total_components = np.sum(Formation_pathway[sim_it][:3])
    Core_frag_frac = Formation_pathway[sim_it][0]/total_components
    Delayed_core_frag_frac = Formation_pathway[sim_it][1]/total_components
    Dynamical_capt_frac = Formation_pathway[sim_it][2]/total_components
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)

#p1 = plt.bar(ind, Dynamical_capt_fracs, 0.95, color='r', hatch='+')
p1 = plt.bar(ind, Core_frag_fracs, 0.95, color='b', linewidth=1, edgecolor='k')#, hatch='+'
p2 = plt.bar(ind, Delayed_core_frag_fracs, 0.95, bottom=Core_frag_fracs, color='m', linewidth=1, edgecolor='k')#, hatch='x'
#p3 = plt.bar(ind, Core_frag_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Dynamical_capt_fracs)), color='b')
p3 = plt.bar(ind, Dynamical_capt_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Core_frag_fracs)), color='r', linewidth=1, edgecolor='k')#, hatch='O'

plt.xlim([-0.6, 5.6])
plt.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=font_size, right=True)
ax.tick_params(axis='both', which='minor', labelsize=font_size, left=True, right=True, top=False, bottom=False)
plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
#plt.yticks(np.arange(0, 1.1,0.1))
ax.tick_params(which='both', direction='in')
plt.xlabel('Initial Gas Mass (M$_\odot$)', fontsize=font_size, labelpad=-0.5)
#ax.set_xticks([1,2,3])
#ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor")
#ax.set_xticklabels(x_labels)#, rotation=45, ha="right", rotation_mode="anchor")
plt.legend((p3[0], p2[0], p1[0]), ('Dynamical capture', 'Delayed core frag.', 'Core fragmentation'), loc='upper right', fontsize=font_size)#, labelspacing=0.2, handletextpad=0.6, borderaxespad=0.3, borderpad=0.2)
plt.ylabel('Fraction', fontsize=font_size, labelpad=-0.5)
plt.ylim([0,1])
plt.savefig('formation_pathway.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

