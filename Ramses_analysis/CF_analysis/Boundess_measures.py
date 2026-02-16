import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import collections

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

pickle_files = ["/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/means_superplot.pkl", "/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/means_superplot.pkl", "/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/means_superplot.pkl", "/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/means_superplot.pkl", "/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/means_superplot.pkl", "/groups/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/means_superplot.pkl"]

for pick_it in range(len(pickle_files)):
    file = open(pickle_files[pick_it], 'rb')
    superplot_dict, Sink_bound_birth, Sink_E_tot, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()
    
    SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE']) - 0.05))
    SFE_5_time = superplot_dict['Times'][SFE_5_ind]
    
    True_new_systems = []
    Mixed_systems = []
    Other_systems = []
    Used_sink_ids = []
    
    #True new systems, It's the first time this sink id appears
    True_core_fragmenation_systems = []
    True_dynamical_capture_systems = []
    True_core_fragmenation_lifetimes = []
    True_dynamical_capture_lifetimes = []
    True_capture_delays = {} #The dealy between the sink formation and when they become bound
    True_fragmentation_delays = {} #This is for true core fragmenation systems, so it looks at the delay between when the components form
    
    #This has some sinks that have previously been found in systems, and sinks that haven't been found in multiple systems
    Mixed_core_fragmentation_systems = []
    Mixed_dynamical_capture_systems = []
    Mixed_core_fragmentation_lifetimes = []
    Mixed_dynamical_capture_lifetimes = []
    Mixed_capture_delays = {} #If the multiple system didn't come into existence as soon as one of the components formed, it is considered a dynamical capture
    Mixed_fragmentation_delays = {} #sometimes a multiple might alread exist but then a new component forms that is bound to the system. This is considered a mixed core fragmenation
    
    #These are systems that don't have any new sink ids. They could have formed from smaller systems becoming capture, or a sink gets ejected from one system and joins another
    Other_core_fragmentation_systems = []
    Other_dynamical_capture_systems = [] #Formed from core fragmenation but the system structure has since changed
    Other_core_fragmentation_lifetimes = []
    Other_dynamical_capture_lifetimes = []
    Other_capture_delays = {}
    Other_fragmentation_delays = {}
    
    for sys_key in Lifetimes_sys.keys():
        sys_formation_time = superplot_dict['System_times'][sys_key][0]
        if sys_formation_time < SFE_5_time:
            sink_ids = flatten(eval(sys_key))
            used_sinks = list(set(sink_ids).intersection(set(Used_sink_ids)))
            if len(used_sinks) == 0:
                True_new_systems.append(sys_key)
                Used_sink_ids = Used_sink_ids + sink_ids
                delay_times = Sink_formation_times[sink_ids].value - superplot_dict['System_times'][sys_key][0]
                if True not in np.array(Sink_bound_birth)[sink_ids]:
                    True_dynamical_capture_systems.append(sys_key)
                    True_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    True_capture_delays.update({sys_key:delay_times})
                elif np.max(delay_times) != 0:
                    True_dynamical_capture_systems.append(sys_key)
                    True_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    True_capture_delays.update({sys_key:delay_times})
                else:
                    True_core_fragmenation_systems.append(sys_key)
                    True_core_fragmenation_lifetimes.append(Lifetimes_sys[sys_key])
                    True_fragmentation_delays.update({sys_key:delay_times})
            elif len(used_sinks) != len(sink_ids):
                Mixed_systems.append(sys_key)
                Used_sink_ids = Used_sink_ids + list(set(used_sinks).symmetric_difference(set(sink_ids)))
                delay_times = Sink_formation_times[sink_ids].value - superplot_dict['System_times'][sys_key][0]
                if True not in np.array(Sink_bound_birth)[sink_ids]:
                    Mixed_dynamical_capture_systems.append(sys_key)
                    Mixed_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    Mixed_capture_delays.update({sys_key:delay_times})
                elif np.max(delay_times) != 0:
                    Mixed_dynamical_capture_systems.append(sys_key)
                    Mixed_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    Mixed_capture_delays.update({sys_key:delay_times})
                else:
                    Mixed_core_fragmentation_systems.append(sys_key)
                    Mixed_core_fragmentation_lifetimes.append(Lifetimes_sys[sys_key])
                    Mixed_fragmentation_delays.update({sys_key:delay_times})
            else:
                Other_systems.append(sys_key)
                delay_times = Sink_formation_times[sink_ids].value - superplot_dict['System_times'][sys_key][0]
                if True not in np.array(Sink_bound_birth)[sink_ids]:
                    Other_dynamical_capture_systems.append(sys_key)
                    Other_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    Other_capture_delays.update({sys_key:delay_times})
                elif np.max(delay_times) != 0:
                    Other_dynamical_capture_systems.append(sys_key)
                    Other_dynamical_capture_lifetimes.append(Lifetimes_sys[sys_key])
                    Other_capture_delays.update({sys_key:delay_times})
                else:
                    Other_dynamical_capture_systems.append(sys_key)
                    Other_core_fragmentation_lifetimes.append(Lifetimes_sys[sys_key])
                    Other_fragmentation_delays.update({sys_key:delay_times})
    import pdb
    pdb.set_trace()
