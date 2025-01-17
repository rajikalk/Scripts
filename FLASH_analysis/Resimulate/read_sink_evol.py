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
            col_tag = row_list
        line_counter = line_counter + 1
        if line_counter > prev_line_counter:
            if row_list[0] != '[00]part_tag':
                if row_list[0] not in sink_data.keys():
                    #If sink id isn't in saved data, add it
                    sink_data.update({row_list[0]:{}})
                    sink_data[row_list[0]].update({col_tag[1].split(']')[-1]:[float(row_list[1])]})
                    sink_data[row_list[0]].update({col_tag[2].split(']')[-1]:[float(row_list[2])]})
                    sink_data[row_list[0]].update({col_tag[3].split(']')[-1]:[float(row_list[3])]})
                    sink_data[row_list[0]].update({col_tag[4].split(']')[-1]:[float(row_list[4])]})
                    sink_data[row_list[0]].update({col_tag[5].split(']')[-1]:[float(row_list[5])]})
                    sink_data[row_list[0]].update({col_tag[6].split(']')[-1]:[float(row_list[6])]})
                    sink_data[row_list[0]].update({col_tag[7].split(']')[-1]:[float(row_list[7])]})
                    sink_data[row_list[0]].update({col_tag[8].split(']')[-1]:[float(row_list[8])]})
                    sink_data[row_list[0]].update({col_tag[9].split(']')[-1]:[float(row_list[9])]})
                    sink_data[row_list[0]].update({col_tag[10].split(']')[-1]:[float(row_list[10])]})
                    sink_data[row_list[0]].update({col_tag[11].split(']')[-1]:[float(row_list[11])]})
                    sink_data[row_list[0]].update({col_tag[12].split(']')[-1]:[float(row_list[12])]})
                    sink_data[row_list[0]].update({col_tag[13].split(']')[-1]:[float(row_list[13])]})
                    sink_data[row_list[0]].update({col_tag[14].split(']')[-1]:[float(row_list[14])]})
                    sink_data[row_list[0]].update({col_tag[15].split(']')[-1]:[float(row_list[15])]})
                else:
                    if float(row_list[1]) in sink_data[row_list[0]][col_tag[1].split(']')[-1]]:
                        #if time in the sink time array
                        match_time = float(row_list[1])
                        remove_keys = []
                        for sink_key in sink_data.keys():
                            #if sink form time is after match time, remove sink
                            if len(sink_data[sink_key]['time']) == 0:
                                remove_keys.append(sink_key)
                            elif sink_data[sink_key]['time'][0] > match_time:
                                remove_keys.append(sink_key)
                        
                        for r_key in remove_keys:
                            del sink_data[r_key]
                            print('removed sink', sink_key)
                                
                        #remove data after this time
                        for sink_key in sink_data.keys():
                            time_ind = np.where(np.array(sink_data[sink_key][col_tag[1].split(']')[-1]]) == match_time)[0][0]
                            for field_key in sink_data[sink_key].keys():
                                sink_data[sink_key][field_key] = sink_data[sink_key][field_key][:time_ind]

                    sink_data[row_list[0]][col_tag[1].split(']')[-1]].append(float(row_list[1]))
                    sink_data[row_list[0]][col_tag[2].split(']')[-1]].append(float(row_list[2]))
                    sink_data[row_list[0]][col_tag[3].split(']')[-1]].append(float(row_list[3]))
                    sink_data[row_list[0]][col_tag[4].split(']')[-1]].append(float(row_list[4]))
                    sink_data[row_list[0]][col_tag[5].split(']')[-1]].append(float(row_list[5]))
                    sink_data[row_list[0]][col_tag[6].split(']')[-1]].append(float(row_list[6]))
                    sink_data[row_list[0]][col_tag[7].split(']')[-1]].append(float(row_list[7]))
                    sink_data[row_list[0]][col_tag[8].split(']')[-1]].append(float(row_list[8]))
                    sink_data[row_list[0]][col_tag[9].split(']')[-1]].append(float(row_list[9]))
                    sink_data[row_list[0]][col_tag[10].split(']')[-1]].append(float(row_list[10]))
                    sink_data[row_list[0]][col_tag[11].split(']')[-1]].append(float(row_list[11]))
                    sink_data[row_list[0]][col_tag[12].split(']')[-1]].append(float(row_list[12]))
                    sink_data[row_list[0]][col_tag[13].split(']')[-1]].append(float(row_list[13]))
                    sink_data[row_list[0]][col_tag[14].split(']')[-1]].append(float(row_list[14]))
                    sink_data[row_list[0]][col_tag[15].split(']')[-1]].append(float(row_list[15]))
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

#plot spin evolution
plt.clf()
for sink_id in write_sink_data.keys():
    plot_time = (write_sink_data[sink_id]['time'] - write_sink_data[sink_id]['time'][0])/3.154e+7
    L_tot = np.sqrt(write_sink_data[sink_id]['anglx']**2 + write_sink_data[sink_id]['angly']**2 + write_sink_data[sink_id]['anglz']**2)
    plt.plot(plot_time, L_tot, label=sink_id)
plt.xlabel('Time (yr)')
plt.ylabel('L (gcm$^2$/s)')
plt.savefig(sink_evol_pickle.split('.pkl')[0] + '.png')
