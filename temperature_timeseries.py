#Read csv

from yt.mods import *
import csv
import matplotlib.pyplot as plt

times = []
temp = []
it = 0

def _Temperature(field, data):
    #This find the temperature, assuming the gas is idea (ie gamma = 5/3)
    #Thermal energy is given in erg/g,
    #Density is given in g/cm^3
    #Number density is given in 1/cm^3
    gamma = 5./3.
    Boltzmann_constant = 1.3806488*(10**(-16)) #erg/K
    top= data["ThermalEnergy"] * data["Density"] * (gamma - 1.)
    bottom = data["NumberDensity"] * Boltzmann_constant
    temperature = top/bottom
    inf = np.all(np.isfinite(data["NumberDensity"]))
    return temperature
    
add_field("Temperature", function=_Temperature, units=r"K")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/highres_hot2/DD00*/CE00*.hierarchy")

#save directory
#save_directory = '/media/DATA/YT_output/hot2/xy-plane-Temperature/'
f = open('temperature_timeseries.csv', 'r+')
f.write('Time, Average_Temperature \n')
it = 0

for pf in ts:
    time = pf.current_time*pf["years"]
    times.append(time)
    
    dd = pf.h.all_data()
    T = dd['Temperature']
    average = sum(T)/len(T)
    temp.append(average)
    
    f.write(str(time) + ',' + str(average) + '\n')
    print time, average
        

plt.clf()
p = plt.plot(times, temp)
plt.xlabel("Time (years)")
plt.ylabel("Average Temperatue ($K$)")
plt.savefig('temperature_timeseries.png')
