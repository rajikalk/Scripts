import pickle
import matplotlib.pyplot as plt
import yt
import numpy as np

pickle_files = ['5_au.pkl', '10_au.pkl']

plt.clf()

for pickle_file in pickle_files:
    file = open(pickle_file, 'rb')
    Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Mass_all, Separation = pickle.load(file)
    file.close()
    
    rad_vel = yt.YTArray(Mean_rad_vel, 'cm/s')
    plt.plot(Time_array, rad_vel.in_units('km/s'), label=pickle_file.split('.')[0])

plt.xlim([0, 10000])
plt.xlabel('Time (yr)')
plt.ylabel('Radial Velocity (km/s)')
plt.legend()
plt.savefig('radial_velocity.png')
