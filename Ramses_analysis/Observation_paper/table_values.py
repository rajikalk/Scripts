import csv
import numpy as np
import sys
import pickle

read_file = sys.argv[1]
pickle_files = []
times = []

with open(read_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0] != '#':
            pickle_file = row[0]
            time = int(row[1])
            pickle_files.append(pickle_file)
            times.append(time)
        

for pit in range(len(pickle_files)):
    file_open = open(pickle_files[pit], 'rb')
    reduced_systems_data = pickle.load(file_open)
    file_open.close()
    
    time_ind = np.argmin(abs(reduced_systems_data['time'].value-times[pit]))
    print("found time", reduced_systems_data['time'][time_ind], "for time", times[pit])
    separation = reduced_systems_data['separation'][0][time_ind]
    inclination = np.rad2deg(np.arcsin(200/separation.value))
    phase = reduced_systems_data['orbital_phase'][0][time_ind]*360
    eccentricity = reduced_systems_data['eccentricity'][0][time_ind]
    period = reduced_systems_data['period'][0][time_ind]
    
    print_values = []
    separation_str = str(int(np.round(separation/10)*10))
    print_values.append(separation_str)
    inclination_str = str(int(np.round(inclination)))
    print_values.append(inclination_str)
    phase_str = str(int(np.round(phase)))
    print_values.append(phase_str)
    eccentricity_str = str(np.round(eccentricity*100)/100)
    print_values.append(eccentricity_str)
    period_str = str(int(np.round(period/100)*100))
    print_values.append(period_str)
    latex_string = ' & '.join(print_values)
    print(latex_string)
