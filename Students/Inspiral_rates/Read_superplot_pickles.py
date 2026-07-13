#!/usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import sys

System_pickle = sys.argv[1]
Sink_particle_pickle = sys.argv[2]
#EDIT TO READ RAW DATA AS WELL
print('Reading pickle')
file_open = open(System_pickle, 'rb')
stuff = pickle.load(file_open)
file_open.close()
print('Finished reading pickle')

Systems_dict = stuff[0]

#plot separations!
print('Plotting system separations')
for key in Systems_dict['System_times'].keys():
    plt.clf()
    plt.semilogy(Systems_dict['System_times'][key], Systems_dict['System_seps'][key])
    plt.xlabel('Time (yr)')
    plt.ylabel('Separation (au)')
    plt.savefig('Separation_vs_time'+key+'.png')

