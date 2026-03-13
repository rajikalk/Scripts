#!/usr/bin/env python

import glob
import yt
import my_ramses_fields_short as myf

file = sorted(glob.glob('data/*/info*txt'))[0]

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}

ds = yt.load(file, units_override=units_override)
dd = ds.all_data()
 
