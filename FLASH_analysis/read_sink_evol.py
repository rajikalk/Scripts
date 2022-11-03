import csv
import sys
import numpy as np
import pickle

sink_evol_file = sys.argv[1]
sink_evol_pickle = sys.argv[2]

sink_data = {}

with open(sink_evol_file, 'r') as f:
    reader = csv.reader(f, delimiter=' ')
    line_counter = 0
    for row in reader:
        line_counter = line_counter + 1
        row_list = [x for x in row if x]
        if row_list[0] == '[00]part_tag':
            col_tag = row_list
        if row_list[0] != '[00]part_tag':
            if row_list[0] not in sink_data.keys():
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
                    match_time = float(row_list[1])
                    for sink_key in sink_data.keys():
                        #if sink form time is after match time, remove sink
                        if sink_data[sink_key]['time'][0] > match_time:
                            del sink_data[sink_key]
                            print('removed sink', sink_key)
                            break
                            
                    #remove data after this time
                    for sink_key in sink_data.keys():
                        if len(np.where(np.array(sink_data[sink_key][col_tag[1].split(']')[-1]]) == match_time)[0]) == 0:
                            import pdb
                            pdb.set_trace()
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
                sink_data[row_list[0]][col_tag[12].split(']')[-1]].append(float(row_list[11]))
                sink_data[row_list[0]][col_tag[13].split(']')[-1]].append(float(row_list[11]))
                sink_data[row_list[0]][col_tag[14].split(']')[-1]].append(float(row_list[11]))
                sink_data[row_list[0]][col_tag[15].split(']')[-1]].append(float(row_list[11]))
            if np.remainder(line_counter, 1000) == 0:
                print('Read up to line', line_counter)
f.close()

for sink_id in sink_data.keys():
    sort_inds = np.argsort(sink_data[sink_id]['time'])
    for para_key in sink_data[sink_id]:
        sink_data[sink_id][para_key] = np.array(sink_data[sink_id][para_key])[sort_inds]

pickle_file = open(sink_evol_pickle, 'wb')
pickle.dump((sink_data), pickle_file)
pickle_file.close()
