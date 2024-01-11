#!/usr/bin/env python
import sys
import glob
import my_flash_module as mym


input_dir = sys.argv[1]
save_dir = sys.argv[2]
plot_time = yt.YTQuantity(4220, 'yr')

files = sorted(glob.glob(input_dir + '*plt_cnt*'))
fn = mym.find_files([plot_time], files)
part_file = 'part'.join(fn.split('plt_cnt'))
ds = yt.load(fn, particle_filename=part_file)

import pdb
pdb.set_trace()

