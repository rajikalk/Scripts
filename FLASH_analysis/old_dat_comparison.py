import csv
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt

Labels = ['Single', 'Mach_0.1-Primary', 'Mach_0.1-Secondary', 'Mach_0.2-Primary', 'Mach_0.2-Secondary']
dat_files = ['/Users/reggie/Documents/Simulation_analysis/FLASH/FLASH_sink_evol/Single_star/sinks_evol.dat', '/Users/reggie/Documents/Simulation_analysis/FLASH/FLASH_sink_evol/Mach_0.1/Lref_10.dat', '/Users/reggie/Documents/Simulation_analysis/FLASH/FLASH_sink_evol/Mach_0.2/Lref_10.dat']
pickle_files = ['Single_mach_0.0.pkl', 'Binary_mach_0.1.pkl', 'Binary_mach_0.2.pkl']

Time_all = []
L_all = []

for dat_it in range(len(dat_files)):

    sink_evol_file = dat_files[dat_it]
    sink_evol_pickle = pickle_files[dat_it]
    sink_data = {}

    with open(sink_evol_file, 'r') as f:
        print("reading ", sink_evol_file)
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
    f.close()

    for sink_id in sink_data.keys():
        sort_inds = np.argsort(sink_data[sink_id]['time'])
        for para_key in sink_data[sink_id]:
            sink_data[sink_id][para_key] = np.array(sink_data[sink_id][para_key])[sort_inds]

    pickle_file = open(sink_evol_pickle, 'wb')
    pickle.dump((sink_data), pickle_file)
    pickle_file.close()

    #plot spin evolution
    plt.clf()
    for sink_id in sink_data.keys():
        L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        time_array = sink_data[sink_id]['time']-sink_data[sink_id]['time'][0]
        Time_all.append(time_array)
        L_all.append(L_tot)
        
        plt.semilogy(sink_data[sink_id]['time'], L_tot, label=sink_id)
        
    plt.xlabel('Time (s)')
    plt.ylabel('L (gcm$^2$/s)')
    plt.savefig(sink_evol_pickle.split('.pkl')[0] + '.png')
    

plt.clf()
for time_it in range(len(Time_all)):
    plt.plot(Time_all[time_it], L_all[time_it], label=Labels[time_it])
plt.xlabel('Time (yr)')
plt.ylabel('L')
plt.legend()
plt.savefig('old_comp.png')
