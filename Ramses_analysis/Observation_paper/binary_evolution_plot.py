import csv
import numpy as np
import sys
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-smooth_acc", "--smooth_accretion", help="do you want to smooth the accrtion", type=str, default='True')
    parser.add_argument("-window", "--smoothing_window", help="how big do you want the smoothing window to be, in term of phase", type=float, default=0.1)
    parser.add_argument("-savename", "--save_image_name", help="What would you like the plot to be saved as?", type=str, default="binary_evolution")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()
read_file = sys.argv[1]
labels = ['B1 (0.8AU)', 'B1', 'B2', 'B3']
linestyles = [':', '-', '-', '-']
colors = ['b', 'b', 'orange', 'g']
if args.smooth_accretion == 'True':
    use_smoothed_mean = True
else:
    use_smoothed_mean = False

if read_file.split('.')[-1] == 'csv':
    pickle_files = []
    times = []

    with open(read_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0][0] != '#':
                pickle_file = row[0]
                time = row[1].split(' ')
                pickle_files.append(pickle_file)
                times.append(time)
        

    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 10.)
    gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0]) #Separation
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #Eccentricity

    Orb_Times_All = []
    Match_Time_All = []
    Match_Sep_All = []
    Match_Acc_All = []
    Match_Ecc_All = []
    Orb_Mean_All = []
    Separation_All = []
    Accretion_All = []
    Eccentricity_All = []

    start_peri_time = 10000

    time_div_number = 10000
    for pit in range(len(pickle_files)):
        file_open = open(pickle_files[pit], 'rb')
        reduced_systems_data = pickle.load(file_open)
        file_open.close()
        
        start_peri_ind = np.argmin(abs(reduced_systems_data['time'].value - start_peri_time))
        time_arr = reduced_systems_data['time']
        sep_arr = reduced_systems_data['separation'][0]
        mdot_arr = reduced_systems_data['mdot'][0]
        ecc_arr = reduced_systems_data['eccentricity'][0]
        
        match_point_x = []
        match_point_x_acc = []
        match_point_y_sep = []
        match_point_y_acc = []
        match_point_y_ecc = []
        for time in times[pit]:
            match_ind = np.argmin(abs(time_arr.value - float(time)))
            match_point_x.append(time_arr[match_ind])
            match_point_y_sep.append(sep_arr[match_ind])
            match_point_x_acc.append(time_arr[match_ind])
            match_point_y_acc.append(mdot_arr[match_ind])
            match_point_y_ecc.append(ecc_arr[match_ind])
        
        left_grad = (reduced_systems_data['separation'][0][start_peri_ind+1:-1]-reduced_systems_data['separation'][0][start_peri_ind:-2])/(reduced_systems_data['time'][start_peri_ind+1:-1]-reduced_systems_data['time'][start_peri_ind:-2])
        right_grad = (reduced_systems_data['separation'][0][start_peri_ind+2:]-reduced_systems_data['separation'][0][start_peri_ind+1:-1])/(reduced_systems_data['time'][start_peri_ind+2:]-reduced_systems_data['time'][start_peri_ind+1:-1])
        if pit == 0:
            extrema_inds = np.argwhere((abs(left_grad*right_grad)>1e-9)&((left_grad*right_grad)<0)) + start_peri_ind
        else:
            extrema_inds = np.argwhere((left_grad*right_grad)<0) + start_peri_ind
        periastron_inds = extrema_inds[::2]
        apastron_inds = extrema_inds[1::2]
        periastron_times = reduced_systems_data['time'][periastron_inds]
        apastron_times = reduced_systems_data['time'][apastron_inds]
        
        orbital_times = []
        peri_it = 0
        ap_it = -1
        match_point_x_orb = []
        for time in reduced_systems_data['time']:
            try:
                if time == apastron_times[ap_it]:
                    peri_it = peri_it + 1
                if time == periastron_times[peri_it]:
                    ap_it = ap_it + 1
                if time >= periastron_times[peri_it] and time < apastron_times[ap_it]:
                    time_scale = apastron_times[ap_it] - periastron_times[peri_it]
                    orb_scale = (time-periastron_times[peri_it])/time_scale/2
                    orb_time = orb_scale + peri_it + 1
                elif time >= apastron_times[ap_it] and time < periastron_times[peri_it]:
                    peri_it = ap_it + 1
                    time_scale = periastron_times[peri_it] - apastron_times[ap_it]
                    orb_scale = (time-apastron_times[ap_it])/time_scale/2
                    orb_time = orb_scale + ap_it + 1.5
                
                if time < periastron_times[peri_it] and peri_it == 0:
                    orbital_times.append(np.nan)
                else:
                    orbital_times.append(orb_time[0].value)
                    
                if time in match_point_x:
                    match_point_x_orb.append(orb_time)
            except:
                pass
                
        window = args.smoothing_window
        mean_orb_time = []
        mean_mdot_arr = []
        for orb_time in orbital_times:
            if np.isnan(orb_time):
                mean_orb_time.append(np.nan)
                mean_mdot_arr.append(np.nan)
            else:
                start_orb = orb_time - window/2
                if start_orb < 1:
                    start_orb = 1
                end_ord = orb_time + window/2
                if end_ord > orbital_times[-1]:
                    end_ord = orbital_times[-1]
                start_ind = np.argmin(abs(np.nan_to_num(orbital_times)-start_orb))
                end_ind = np.argmin(abs(np.nan_to_num(orbital_times)-end_ord))
                mean_mdot = np.mean(mdot_arr[start_ind:end_ind])
                mean_orb = np.mean(orbital_times[start_ind:end_ind])
                mean_orb_time.append(mean_orb)
                mean_mdot_arr.append(mean_mdot)
                
        if use_smoothed_mean == True:
            match_point_x_acc = []
            match_point_y_acc = []
            for time in times[pit]:
                match_ind = np.argmin(abs(time_arr.value - float(time)))
                match_point_x_acc.append(mean_orb_time[match_ind])
                match_point_y_acc.append(mean_mdot_arr[match_ind])
        
        '''
        mdot_mean_arr = []
        time_mean_arr = []
        for time_it in range(len(time_arr)):
            start_time = time_arr[time_it].value-window/2
            if start_time < time_arr[0].value:
                start_time = time_arr[0].value
            end_time = time_arr[time_it].value+window/2
            if end_time > time_arr[-1].value:
                end_time = time_arr[-1].value
            start_ind = np.argmin(abs(time_arr.value - start_time))
            end_ind = np.argmin(abs(time_arr.value - end_time))
            mean_mdot = np.mean(mdot_arr[start_ind:end_ind])
            mean_time = np.mean(time_arr[start_ind:end_ind])
            mdot_mean_arr.append(mean_mdot)
            time_mean_arr.append(mean_time/time_div_number)
        '''
            
        time_arr = reduced_systems_data['time']/time_div_number
        
        #trim_arrays
        start_trim_ind = np.argmin(abs(np.nan_to_num(orbital_times)-1))
        end_trim_ind = np.argmin(abs(np.nan_to_num(orbital_times)-10))
        orbital_times = orbital_times[start_trim_ind:end_trim_ind+1]
        sep_arr = sep_arr[start_trim_ind:end_trim_ind+1]
        mdot_arr = mdot_arr[start_trim_ind:end_trim_ind+1]
        mean_orb_time = mean_orb_time[start_trim_ind:end_trim_ind+1]
        mean_mdot_arr = mean_mdot_arr[start_trim_ind:end_trim_ind+1]
        ecc_arr = ecc_arr[start_trim_ind:end_trim_ind+1]
            
        
        #ax1.semilogy(time_arr, sep_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        ax1.semilogy(orbital_times, sep_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        ax1.scatter(match_point_x_orb, match_point_y_sep, marker='o', color=colors[pit], facecolors='none')
        
        #ax2.semilogy(time_mean_arr, mdot_mean_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        if use_smoothed_mean == True:
            ax2.semilogy(mean_orb_time, np.array(mean_mdot_arr)*100000, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        else:
            ax2.semilogy(orbital_times, mdot_arr*100000, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        ax2.scatter(match_point_x_orb, np.array(match_point_y_acc)*100000, marker='o', color=colors[pit], facecolors='none')
        
        #ax3.semilogy(time_arr, ecc_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        #ax3.semilogy(orbital_times, ecc_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        ax3.plot(orbital_times, ecc_arr, label=labels[pit], ls=linestyles[pit], linewidth=0.5, alpha=0.5, color=colors[pit])
        ax3.scatter(match_point_x_orb, match_point_y_ecc, marker='o', color=colors[pit], facecolors='none')
        
        Orb_Times_All.append(orbital_times)
        Separation_All.append(sep_arr)
        Match_Time_All.append(match_point_x_orb)
        Match_Sep_All.append(match_point_y_sep)
        Match_Acc_All.append(match_point_y_acc)
        Match_Ecc_All.append(match_point_y_ecc)
        Orb_Mean_All.append(mean_orb_time)
        Accretion_All.append(mean_mdot_arr)
        Eccentricity_All.append(ecc_arr)
    
    file_open = open(args.save_image_name+'.pkl', 'wb')
    pickle.dump((Orb_Times_All, Separation_All, Match_Time_All, Match_Sep_All, Match_Acc_All, Match_Ecc_All, Orb_Mean_All, Accretion_All, Eccentricity_All),file_open)
    file_open.close()
        
elif read_file.split('.')[-1] == 'pkl':
    file_open = open(read_file, 'rb')
    Orb_Times_All, Separation_All, Match_Time_All, Match_Sep_All, Match_Acc_All, Match_Ecc_All, Orb_Mean_All, Accretion_All, Eccentricity_All = pickle.load(file_open)
    file_open.close()
    
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 10.)
    gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0]) #Separation
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #Accretion Rate
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #Eccentricity
    
    for it in range(len(Orb_Times_All)):
        ax1.semilogy(Orb_Times_All[it], Separation_All[it], label=labels[it], ls=linestyles[it], linewidth=0.75, alpha=0.75, color=colors[it])
        
        ax2.semilogy(Orb_Mean_All[it], np.array(Accretion_All[it])*100000, label=labels[it], ls=linestyles[it], linewidth=0.75, alpha=0.75, color=colors[it])
        
        ax3.plot(Orb_Times_All[it], Eccentricity_All[it], label=labels[it], ls=linestyles[it], linewidth=0.75, alpha=0.75, color=colors[it])
    
    for it in range(len(Orb_Times_All)):
        ax1.scatter(Match_Time_All[it], Match_Sep_All[it], marker='o', color=colors[it], facecolors='none')
        ax2.scatter(Match_Time_All[it], np.array(Match_Acc_All[it])*100000, marker='o', color=colors[it], facecolors='none')
        ax3.scatter(Match_Time_All[it], Match_Ecc_All[it], marker='o', color=colors[it], facecolors='none')

ax1.set_ylabel('Separation (AU)')
plt.setp([ax1.get_xticklabels() for ax1 in fig.axes[:-1]], visible=False)
#ax1.set_ylim([50, 3e3])
ax2.set_ylabel('Accretion ($10^{-5}$M$_\odot$/yr)')
plt.setp([ax2.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
#ax2.set_ylim([2.5e-6, 6e-5])
ax3.set_ylabel('Eccentricity')
ax3.set_xlabel('Orbit #')
#ax3.set_ylim([0.3,1.5])
ax3.set_xticks(np.arange(12))
ax1.legend(loc='best')
ax1.set_xlim([1, 10])

ax1.xaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax1.grid(which='both', alpha=0.2)
ax1.grid(which='minor', alpha=0.2, linestyle='--')
ax2.grid(which='both', alpha=0.2)
ax2.grid(which='minor', alpha=0.2, linestyle='--')
ax2.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
ax2.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax3.grid(which='both', alpha=0.2)
ax3.grid(which='minor', alpha=0.2, linestyle='--')
ax3.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
ax3.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.tick_params(axis='y', which='major', direction="in")
ax1.tick_params(axis='y', which='minor', direction="in")
ax1.tick_params(axis='x', which='major', direction="in")
ax1.tick_params(axis='x', which='minor', direction="in")
ax2.tick_params(axis='y', which='major', direction="in")
ax2.tick_params(axis='y', which='minor', direction="in")
ax2.tick_params(axis='x', which='major', direction="in")
ax2.tick_params(axis='x', which='minor', direction="in")
ax3.tick_params(axis='y', which='major', direction="in")
ax3.tick_params(axis='y', which='minor', direction="in")
ax3.tick_params(axis='x', which='major', direction="in")
ax3.tick_params(axis='x', which='minor', direction="in")

plt.savefig(args.save_image_name+'.pdf', format='pdf', bbox_inches='tight')
    
    
