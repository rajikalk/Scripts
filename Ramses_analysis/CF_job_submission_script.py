import numpy as np
import os
import subprocess

dirs = ["/groups/astro/rlk/Analysis_plots/Ramses/Global/G400", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G200", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/256", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G50"]
Match_Methods = ["Ratio_vis_to_all", "Number_Stars_50", "Number_Stars", "SFE", "SFE_t_ff", "M_tot_150", "Ratio_vis_to_all"]#, "Number_Stars_50"]
Boundness = ["Unbound", "Bound"]
L_Limits = ["L_Limits", "No_Limits"]
global_pickles = ["/groups/astro/rlk/Analysis_plots/Ramses/Global/G400/stars_imf_G400.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G200/stars_imf_G200.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/256/stars_imf_G100.pkl", "/groups/astro/rlk/Analysis_plots/Ramses/Global/G50/stars_imf_G50.pkl"]
n_nodes = [6, 3, 3, 2, 1]

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-rm_pkl", "--remove_pickle", help="do ou want to remove existing pickles and star afresh?", type=str, default='True')
    parser.add_argument("-use_t_s", "--use_t_spread", help="do you want to use t_spread?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
args = parse_inputs()
nit = -1
for dir in dirs:
    nit = nit + 1
    for method in Match_Methods:
        for bound_state in Boundness:
            for limits in L_Limits:
                save_dir = dir + '/' + method + '/' + bound_state + '/' + limits + '/'
                #if method == "Number_Stars":
                #    save_dir = save_dir + 'No_Obs_Limits/'
                if os.path.exists(save_dir) == False:
                    os.makedirs(save_dir)
                os.chdir(save_dir)
                if args.remove_pickle == 'True':
                    subprocess.call('rm *.pkl', shell=True)
                    subprocess.call('rm slurm*.out', shell=True)
                    subprocess.call('rm *.jpg', shell=True)
                    subprocess.call('rm *.pdf', shell=True)
                    subprocess.call('rm *.png', shell=True)
                if method == 'SFE':
                    job_string = 'SFE'
                    method_ind = 1
                elif method == 'SFE_t_ff':
                    job_string = 'STF'
                    method_ind = 2
                elif method == 'M_tot_150':
                    job_string = 'MTo'
                    method_ind = 4
                elif method == 'Ratio_vis_to_all':
                    job_string = 'Rat'
                    method_ind = 5
                else:
                    method_ind = 3
                    if '_50' in method:
                        n_vis_thres = 50
                        job_string = 'NS50'
                    else:
                        n_vis_thres = 115
                        job_string = 'NoS'
                if 'G100' in dir:
                    job_tag = limits[0] + bound_state[0] + job_string + dir.split('Global/')[1][:2] + dir.split('G100/')[1]
                else:
                    job_tag = limits[0] + bound_state[0] + job_string + dir.split('Global/')[1][:3]
                nodes = n_nodes[nit]
    
                use_short = bool(np.random.randint(0,2))
                add_2 = bool(np.random.randint(0,2))
                queue = 'astro'
                t_per_n = 20
                if add_2:
                    queue = queue + '2'
                    t_per_n = 48
                if '2' not in queue:
                    nodes = nodes*2
                if use_short:
                    queue = queue + '_short'
                else:
                    queue = queue + '_long'
                time_str = '6:00:00'
                '''
                if queue == 'astro2_short':
                    time_str = '12:00:00'
                elif queue == 'astro_short':
                    time_str = '1-00:00:00'
                else:
                    time_str = '5-00:00:00'
                '''
                    
                if bound_state == 'Bound':
                    Bound_bool = "True"
                else:
                    Bound_bool = "False"
                
                run_string = "srun python /groups/astro/rlk/Scripts/Ramses_analysis/multiplicity_fraction_res_per_bin_13.py /lustre/astro/troels/IMF_512/binary_analysis/data/ ./ -pickle cf_hist -verbose True -global_data "+global_pickles[nit]+" -bound "+str(Bound_bool) + " -match_meth " + str(method_ind)
                
                if args.use_t_spread != 'True':
                    run_string = run_string + "-use_t_s False"
                
                if limits == 'No_Limits':
                    run_string = run_string + " -upper_L 100000"
                
                if method_ind == 3:
                    run_string = run_string + " -use_t_s False -n_vis_thres "+str(n_vis_thres)# -acc_lim 1.e-12 -lower_L 0"
                
                run_string = run_string + " 1>cf.out00 2>&1\n"
        
                f = open('cf_pickles.sh', 'w')
                f.write("#!/bin/bash\n")
                f.write("# SLURM resource specifications\n")
                f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
                f.write("#SBATCH --job-name="+job_tag+"   # shows up in the output of 'squeue'\n")
                f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
                f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
                f.write("#SBATCH --nodes="+str(nodes)+"              # number of nodes allocated for this job\n")
                f.write("#SBATCH --ntasks-per-node="+str(t_per_n)+"    # number of MPI ranks per node\n")
                f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
                f.write("#SBATCH --mail-type=ALL\n")
                
                f.write(run_string)
                f.close()

                subprocess.call('sbatch cf_pickles.sh', shell=True)
                print("submitted job", save_dir+'cf_pickles.sh')
                
                
