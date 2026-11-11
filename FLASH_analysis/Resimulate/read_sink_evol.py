import csv
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
from copy import deepcopy

sink_evol_file = sys.argv[1]
sink_evol_pickle = sys.argv[2]

try:
    if os.path.exists(sink_evol_pickle):
        file = open(sink_evol_pickle, 'rb')
        sink_data, prev_line_counter = pickle.load(file)
        file.close()
        for sink_id in sink_data.keys():
            for para_key in sink_data[sink_id].keys():
                sink_data[sink_id][para_key] = list(sink_data[sink_id][para_key])
    else:
        prev_line_counter = 0
        sink_data = {}
except:
    prev_line_counter = 0
    sink_data = {}

with open(sink_evol_file, 'r') as f:
    print("reading ", sink_evol_file)
    reader = csv.reader(f, delimiter=' ')
    line_counter = 0
    for row in reader:
        row_list = [x for x in row if x]
        if row_list[0] == '[00]part_tag':
            col_tags = row_list
        line_counter = line_counter + 1
        if line_counter > prev_line_counter:
            if row_list[0] != '[00]part_tag':
                if row_list[0] not in sink_data.keys():
                    #If sink id isn't in saved data, add it
                    sink_data.update({row_list[0]:{}})
                    for col_tag in col_tags[1:]:
                        sink_data[row_list[0]].update({col_tag.split(']')[-1]:np.array([float(row_list[col_tags.index(col_tag)])])})
                else:
                    if float(row_list[1]) in sink_data[row_list[0]][col_tags[1].split(']')[-1]]:
                        #Get index of the time that is already in the sink particle dictionary
                        replace_ind = np.where(sink_data[row_list[0]][col_tags[1].split(']')[-1]] == float(row_list[1]))[0][0]
                        
                        #remove that offending index
                        for col_tag in col_tags[1:]:
                            sink_data[row_list[0]][col_tag.split(']')[-1]] = np.delete(sink_data[row_list[0]][col_tag.split(']')[-1]], replace_ind)
                                
                    #Add new data for existing arrays
                    for col_tag in col_tags[1:]:
                        sink_data[row_list[0]][col_tag.split(']')[-1]] = np.append(float(row_list[col_tags.index(col_tag)]), sink_data[row_list[0]][col_tag.split(']')[-1]])

                if np.remainder(line_counter, 5000) == 0:
                    print('Read up to line', line_counter)
                    write_sink_data = deepcopy(sink_data)
                    for sink_id in write_sink_data.keys():
                        sort_inds = np.argsort(write_sink_data[sink_id]['time'])
                        for para_key in write_sink_data[sink_id].keys():
                            write_sink_data[sink_id][para_key] = np.array(write_sink_data[sink_id][para_key])[sort_inds]
                    pickle_file = open(sink_evol_pickle, 'wb')
                    pickle.dump((write_sink_data, line_counter), pickle_file)
                    pickle_file.close()
f.close()

write_sink_data = deepcopy(sink_data)
for sink_id in write_sink_data.keys():
    sort_inds = np.argsort(write_sink_data[sink_id]['time'])
    for para_key in write_sink_data[sink_id].keys():
        write_sink_data[sink_id][para_key] = np.array(write_sink_data[sink_id][para_key])[sort_inds]

pickle_file = open(sink_evol_pickle, 'wb')
pickle.dump((write_sink_data, line_counter), pickle_file)
pickle_file.close()

#plot mass and specific angular momentum evolution
Msun = 1.989e33   # g
yr = 3.154e7      # s

fig, (ax_mass, ax_spin) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))

for sink_id in write_sink_data.keys():
    sink = write_sink_data[sink_id]
    plot_time = (sink['time'] - sink['time'][0])/yr
    L_tot = np.sqrt(sink['anglx']**2 + sink['angly']**2 + sink['anglz']**2)
    L_spec = L_tot/sink['mass']

    ax_mass.plot(plot_time, sink['mass']/Msun, label=sink_id)
    ax_spin.plot(plot_time, L_spec, label=sink_id)

ax_mass.set_ylabel('Mass (M$_\\odot$)')
ax_spin.set_ylabel('Specific angular momentum (cm$^2$/s)')
ax_spin.set_xlabel('Time (yr)')
ax_mass.legend(loc='best', fontsize='small')

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig(sink_evol_pickle.split('.pkl')[0] + '.png', dpi=300)
