import numpy as np
import os
import subprocess

dirs = ["/groups/astro/rlk/Analysis_plots/Ramses/Global/G400", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G200", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/256", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G50"]
Match_Methods = ["Number_Stars", "SFE", "SFE_t_ff", "M_tot_150"]
Boundness = ["Unbound", "Bound"]
L_Limits = ["L_Limits", "No_Limits"]

nit = -1
for dir in dirs:
    nit = nit + 1
    for method in Match_Methods:
        for bound_state in Boundness:
            for limits in L_Limits:
                astro_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/L_mean_hist.jpg'
                save_dir = dir.split('Global/')[1] + '/L_plots/L_hist_' + method + '_' + bound_state + '_' + limits + '.jpg'

                copy_string = 'rsync -vaz astro07-travel:'+astro_dir + ' ' + save_dir
                subprocess.call(copy_string, shell=True)
                print("copied to", save_dir)
                
                if method == "Number_Stars":
                    astro_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/visible_stars.png'
                    save_dir = dir.split('Global/')[1] + '/L_plots/visible_stars_' + bound_state + '_' + limits + '.jpg'
                    copy_string = 'rsync -vaz astro07-travel:'+astro_dir + ' ' + save_dir
                    subprocess.call(copy_string, shell=True)
                    print("copied to", save_dir)
                elif method == "SFE":
                    astro_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/SFE.png'
                    save_dir = dir.split('Global/')[1] + '/L_plots/SFE_' + bound_state + '_' + limits + '.jpg'
                    copy_string = 'rsync -vaz astro07-travel:'+astro_dir + ' ' + save_dir
                    subprocess.call(copy_string, shell=True)
                    print("copied to", save_dir)
                elif method == "SFE_t_ff":
                    astro_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/SFE_n.png'
                    save_dir = dir.split('Global/')[1] + '/L_plots/SFE_n_' + bound_state + '_' + limits + '.jpg'
                    copy_string = 'rsync -vaz astro07-travel:'+astro_dir + ' ' + save_dir
                    subprocess.call(copy_string, shell=True)
                    print("copied to", save_dir)
                else:
                    astro_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/M_total.png'
                    save_dir = dir.split('Global/')[1] + '/L_plots/M_total_' + bound_state + '_' + limits + '.jpg'
                    copy_string = 'rsync -vaz astro07-travel:'+astro_dir + ' ' + save_dir
                    subprocess.call(copy_string, shell=True)
                    print("copied to", save_dir)
                
                    
