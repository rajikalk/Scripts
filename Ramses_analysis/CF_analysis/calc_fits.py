import os
import subprocess

pickle_dir_prefix = '/lustre/astro/rlk/Analysis_plots/Ramses/Superplots/'
pickle_dir_suffix = 'Replace_with_most_massive_star/CF_plots/CF_movies/'
save_dir_root = '/groups/astro/rlk/Analysis_plots/Fit_over_time_plots/'

sub_dirs = ['G50/', 'G100/256/', 'G100/512/', 'G125/', 'G150/', 'G200/', 'G400/']
bound_dirs = ['Bound_all_stars/', 'Unbound_visible_tobin_stars/']
proj_dirs = ['X_proj/', 'Y_proj/', 'Z_proj/']

for sub_dir in sub_dirs:
    for bound_dir in bound_dirs:
        pickle_dir = pickle_dir_prefix + sub_dir + pickle_dir_suffix + bound_dir+'cf_hist'
        save_dir = save_dir_root + sub_dir + bound_dir
        call_string = 'python ~/Scripts/Ramses_analysis/CF_analysis/fit_over_time.py -pickle '+pickle_dir+' -savedir ' + save_dir + ' -t_spread ' + str(100000)
        subprocess.call(call_string, shell=True)
        for proj_dir in proj_dirs:
            pickle_dir = pickle_dir_prefix + sub_dir + pickle_dir_suffix + bound_dir+proj_dir+'cf_hist'
            save_dir = save_dir_root + sub_dir + bound_dir+proj_dir
            call_string = 'python ~/Scripts/Ramses_analysis/CF_analysis/fit_over_time.py -pickle '+pickle_dir+' -savedir ' + save_dir  + ' -t_spread ' + str(100000)
            subprocess.call(call_string, shell=True)
