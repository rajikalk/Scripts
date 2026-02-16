import pickle
import matplotlib.pyplot as plt
import numpy as np

pickle_dir = '/lustre/hpc/astro/rlk/Analysis_plots/Ramses/Obs_comp/'
pickles = ['Sink_49/particle_data_neat.pkl', 'Sink_91/particle_data_neat.pkl', 'Sink_91/High_resolution/particle_data_neat.pkl', 'Sink_164/particle_data_neat.pkl']

highlight_times = [[104300, 98600, 97300], [52100, 50300], [45000], [50400]]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(10, 8), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0.0)

for fit in range(len(pickles)):
    file_open = open(pickles[fit], 'rb')
    particle_data, something, sink_ind = pickle.load(file_open)
    file_open.close()
    sep_x = particle_data['posx'][0]-particle_data['posx'][1]
    sep_y = particle_data['posy'][0]-particle_data['posy'][1]
    sep_z = particle_data['posz'][0]-particle_data['posz'][1]
    sep = np.sqrt(sep_x**2 + sep_y**2 + sep_z**2)
    axs[fit].semilogy(particle_data['time'], sep, alpha=0.5, color='k')
    print('plotted sink', sink_ind)
    for pit in highlight_times[fit]:
        time_it = np.argmin(abs(particle_data['time'].value-pit))
        axs[fit].plot(pit, sep[time_it], marker='o', color='b')
    axs[fit].set_ylabel('Separation (AU)')
plt.xlabel('Time (yr)')
plt.savefig('highlight_times.png')
