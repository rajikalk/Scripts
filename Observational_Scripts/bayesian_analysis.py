import binary_orbit as bo
import numpy as np
from astropy.time import Time

my_orbs = bo.random_orbits()

Object = []
Region = []
Coords = []
Membership_probability = []
IR_excess = []
Disk_visual = []
Magnitudes = []
No_obs = []
SB1_flag = []
SB2_flag = []
RV_variation = []
Temp_sptype = []
Pref_template = []
Obs_info = []
SNR = []

RV_standard_info = {}

#Read in RV standard list
header = 0
with open('/Users/rajikak/Observational_Data/RV_standard_list.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            RV_standard_info[row[0]] = (float(6), float(7))
        else:
            header = 1

print "Reading in current spreadsheet"
header = 0
with open(args.input_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if 'U4' in row[0]:
                row[0] = 'UCAC4' + row[0].split('U4')[-1]
            Object.append(row[0])
            Region.append(row[1])
            Coords.append((row[2], row[3]))
            Membership_probability.append(int(row[4]))
            IR_excess.append(row[5])
            Disk_visual.append(row[6])
            Magnitudes.append([float(row[7]), float(row[8]), float(row[9])])
            No_obs.append(int(row[10]))
            if row[11] == 'FALSE' or row[11] == 'False':
                SB1_flag.append(False)
            else:
                SB1_flag.append(True)
            if row[12] == 'FALSE' or row[12] == 'False':
                SB2_flag.append(False)
            else:
                SB2_flag.append(True)
            RV_variation.append(row[13])
            Pref_template.append(row[14])
            Temp_sptype.append(row[15])
            if len(row) > 16:
                Obs = np.array(row[16:])
                Obs = np.delete(Obs, np.where(Obs==''))
                #if len(Obs) > 5:
                #    Obs = np.reshape(Obs, (len(Obs)/5, 5))
                Obs = np.reshape(Obs, (len(Obs)/6, 6))
                for ind_obs in Obs:
                    if '/' in ind_obs[0]:
                        new_format = '20' + ind_obs[0].split('/')[-1] + '-' + ind_obs[0].split('/')[-2] + '-' + ("%02d" % int(ind_obs[0].split('/')[-3]))
                        ind_obs[0] = new_format
            else:
                Obs = np.array([])
            Obs_info.append(Obs)
        if header == 0:
            header = 1

for obj in range(len(Object)):
    #Get M_1 from the template spectral type
    Pref_template_name = Pref_template[obj].split('_')[0]

