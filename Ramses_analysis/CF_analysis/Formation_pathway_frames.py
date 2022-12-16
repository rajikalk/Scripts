#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
import my_ramses_fields as myf
import csv

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-22')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-19)
    parser.add_argument("-G_mass", "--Gas_mass", type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = '100'

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
mym.set_units(units_override)

#------------------------------
Sim_path = '/lustre/astro/troels/IMF_256_fixed_dt/data/'
files = sorted(glob.glob(Sim_path+"*/info*.txt"))
txt_files = sorted(glob.glob(Sim_path+"*/stars_output.snktxt"))
sim_file_times = []

for output_txt in txt_files:
    with open(output_txt, 'rU') as txt_file:
        reader = csv.reader(txt_file)
        for row in reader:
            time_val = float(row[0].split('   ')[-2])
            sim_file_times.append(time_val)

Interested_sinks = [36, 14, 2]
Other_sink = [4, [10, [5, 9]], [1, 3]]

#Bound core fragmentation pathway
Primary_form_time = 1.0365265956563827
Secondary_form_time = 1.0460617956407776

m_times = [Secondary_form_time, Primary_form_time]

import pdb
pdb.set_trace()


#Delay core frag pathway
Primary_form_time = 1.0387929956526736
Secondary_form_time = 1.040190745650386
System_bound_time = 1.0405948456497247

m_times = [Primary_form_time, Secondary_form_time, System_bound_time]

#Dynamical capture
Star_form_time = 1.0358556956574807
Capture_time = 1.036077995657117

m_times = [Star_form_time, Capture_time]
