import pickle
import sys

sim_id = sys.argv[1]

full_search_pickle = '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/'+sim_id+'/Full_sink_data/Fast_analysis/Full_first_sys_check/sink_birth_all_delayed_core_frag_cleaned.pkl'
first_pickle = '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/'+sim_id+'/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl'

file_open = open(full_search_pickle, 'rb')
Sink_birth_all_full = pickle.load(file_open)
file_open.close()

file_open = open(first_pickle, 'rb')
Sink_birth_all_first = pickle.load(file_open)
file_open.close()

for key in Sink_birth_all_full.keys():
    if Sink_birth_all_full[key] != Sink_birth_all_first[key]:
        import pdb
        pdb.set_trace()


