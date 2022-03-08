import numpy as np
import matplotlib.pyplot as plt

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

#Formation_pathway = [[252, 88, 105], [624, 859, 1209], [1893, 2000, 5105], [1255, 4381, 11458], [2921, 5951, 22602], [2172, 5412, 32412]]
Formation_pathway = [[96, 129, 108], [624, 859, 1209], [1536, 994, 2917], [1340, 1746, 5248], [3353, 2539, 9586], [6559, 4865, 22248]]
Formation_pathway = [[18, 14, 14, 287], [28, 35, 43, 2586], [55, 33, 111, 5248], [54, 56, 135, 8089], [104, 71, 242, 15061], [np.nan, np.nan, np.nan, np.nan]]

x_labels = ['Core frag.', 'Delayed core frag.', 'Dynamical capture']

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
Core_frag_fracs = []
Delayed_core_frag_fracs = []
Dynamical_capt_fracs = []
Other_fracs = []
ind = np.arange(len(subplot_titles))

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
plt.xlabel('Initial Gas Mass (M$_\odot$)')
#ax.set_xticks([1,2,3])
#ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor")
#ax.set_xticklabels(x_labels)#, rotation=45, ha="right", rotation_mode="anchor")
plt.legend((p4[0], p3[0], p2[0], p1[0]), ('Core fragmentation', 'Delayed core fragmentation', 'Dynamical capture', 'Other'), loc='lower right')
plt.ylabel('Fraction')
plt.ylim([0,1])
plt.savefig('formation_pathway_other.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
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

p1 = plt.bar(ind, Dynamical_capt_fracs, 0.95, color='r')
p2 = plt.bar(ind, Delayed_core_frag_fracs, 0.95, bottom=Dynamical_capt_fracs, color='m')
p3 = plt.bar(ind, Core_frag_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Dynamical_capt_fracs)), color='b')

plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
#plt.yticks(np.arange(0, 1.1,0.1))
ax.tick_params(which='both', direction='in')
plt.xlabel('Initial Gas Mass (M$_\odot$)')
#ax.set_xticks([1,2,3])
#ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor")
#ax.set_xticklabels(x_labels)#, rotation=45, ha="right", rotation_mode="anchor")
plt.legend((p3[0], p2[0], p1[0]), ('Core fragmentation', 'Delayed core fragmentation', 'Dynamical capture'), loc='lower right')
plt.ylabel('Fraction')
plt.ylim([0,1])
plt.savefig('formation_pathway.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

