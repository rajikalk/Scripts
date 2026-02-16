import os
import numpy as np
import subprocess

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-make_cf_his", "--make_cf_histogram", help="do you need to make the cf_hist.pkl?", type=str, default='True')
    parser.add_argument("-make_movie", "--make_movie_frames", help="do you need to make the movie frames", type=str, default='False')
    parser.add_argument("-rm_frames", "--remove_existing_frames", help="Do you want replot the frames", type=str, default='False')
    parser.add_argument("-rm_pickle", "--remove_existing_pickle", help="Do you want remake the pickle", type=str, default='False')
    parser.add_argument("-short_only", "--short_queues_only", help="Do you want to only use the short queues?", type=str, default='False')
    parser.add_argument("-astro02", "--astro02_only", help="Do you want to only use the astro02 queue?", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args


args = parse_inputs()
root_dir = '/lustre/hpc/astro/rlk/rlk/Analysis_plots/Ramses/Superplots/'
top_dirs = ['G50/', 'G100/256/', 'G100/512/', 'G125/', 'G150/', 'G200/', 'G400/']
mid_dirs = '2021_Tobin_analysis_update/'
bot_dirs = ['L_limits/', 'No_L_limits/']
proj_dirs = ['X_proj/', 'Y_proj/', 'Z_proj/']
global_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G50/stars_imf_G50.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G100/256/stars_imf_G100.pkl",
    "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl",
    "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G125/stars_imf_G125.pkl",
    "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G150/stars_imf_G150.pkl",
    "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G200/stars_imf_G200.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Ramses/Global/G400/stars_imf_G400.pkl"]
job_tag = ['G50', 'G1256', 'G1512', 'G125', 'G150', 'G200', 'G400']
n_nodes = [3, 4, 5, 5, 10, 5, 5, 10]
#n_nodes = [1, 2, 3, 3, 5, 3, 3, 6]

job_it = -1
for tdir in top_dirs:
    job_it = job_it + 1
    for bdir in bot_dirs:
        save_dir = root_dir+tdir + mid_dirs + bdir
        if args.make_cf_histogram == 'True':
            if os.path.exists(save_dir) == False:
                os.makedirs(save_dir)
            os.chdir(save_dir)
            if args.remove_existing_pickle != 'False':
                subprocess.call('rm *.pkl', shell=True)
            job_string = bdir[0] + job_tag[job_it]
            nodes = n_nodes[job_it]
            if args.short_queues_only != 'False':
                use_short = True
            else:
                use_short = bool(np.random.randint(0,2))
            if args.astro02_only != 'False':
                add_2 = True
            else:
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
            if queue == 'astro2_short':
                time_str = '12:00:00'
            elif queue == 'astro_short':
                time_str = '1-00:00:00'
            elif queue == 'astro2_long':
                time_str = '5-00:00:00'
            else:
                time_str = '10-00:00:00'

            if bdir.split('_')[0] == 'Bound':
                L_bool = "True"
            else:
                L_bool = "False"
                
            run_string = "srun python /groups/astro/rlk/Scripts/Ramses_analysis/CF_analysis/multiplicity_fraction_res_per_bin_13.py /lustre/astro/troels/IMF_512/binary_analysis/data/ ./ -pickle cf_hist -verbose True -global_data " + global_pickles[job_it] + " -bound True -entire_sim True "
            if L_bool == "False":
                run_string = run_string + "-lower_L 0 -upper_L 10000 -acc_lim 0 -replace_ind_type M "
                
            run_string = run_string + "1>cf.out00 2>&1\n"
            
            f = open('cf_pickles.sh', 'w')
            f.write("#!/bin/bash\n")
            f.write("# SLURM resource specifications\n")
            f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
            f.write("#SBATCH --job-name="+job_string+"   # shows up in the output of 'squeue'\n")
            f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
            f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
            f.write("#SBATCH --nodes="+str(nodes)+"              # number of nodes allocated for this job\n")
            f.write("#SBATCH --ntasks-per-node="+str(t_per_n)+"    # number of MPI ranks per node\n")
            f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
            f.write("#SBATCH --mail-type=ALL\n")
            
            f.write(run_string)
            f.close()

            if os.path.exists(save_dir+"cf_hist.pkl") == False:
                subprocess.call('sbatch cf_pickles.sh', shell=True)
                print("submitted job", save_dir+'cf_pickles.sh')
        
        if args.make_movie_frames == 'True' and os.path.exists(save_dir+'cf_pickles.sh'):
            if os.path.exists(save_dir) == False:
                os.makedirs(save_dir)
            os.chdir(save_dir)
            if args.remove_existing_frames != 'False':
                subprocess.call('rm movie_frame*.jpg', shell=True)
            job_string = "M" + bdir[0] + job_tag[job_it]
            nodes = n_nodes[-1*job_it]
            if args.short_queues_only != 'False':
                use_short = True
            else:
                use_short = bool(np.random.randint(0,2))
            if args.astro02_only != 'False':
                add_2 = True
            else:
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
            if queue == 'astro2_short':
                time_str = '12:00:00'
            elif queue == 'astro_short':
                time_str = '1-00:00:00'
            elif queue == 'astro2_long':
                time_str = '5-00:00:00'
            else:
                time_str = '10-00:00:00'
                
            run_string = "srun python ~/Scripts/Ramses_analysis/CF_analysis/cf_movie.py ./ ./ -pickle cf_hist.pkl 1>cf_movie.out00 2>&1\n"
            
            f = open('cf_movie.sh', 'w')
            f.write("#!/bin/bash\n")
            f.write("# SLURM resource specifications\n")
            f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
            f.write("#SBATCH --job-name="+job_string+"   # shows up in the output of 'squeue'\n")
            f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
            f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
            f.write("#SBATCH --nodes="+str(nodes*2)+"              # number of nodes allocated for this job\n")
            f.write("#SBATCH --ntasks-per-node="+str(t_per_n)+"    # number of MPI ranks per node\n")
            f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
            f.write("#SBATCH --mail-type=ALL\n")
            
            f.write(run_string)
            f.close()

            if os.path.exists(save_dir+"cf_hist.pkl"):# and os.path.exists(save_dir+"cf_movie.mp4")==False:
                subprocess.call('sbatch cf_movie.sh', shell=True)
                print("submitted job", save_dir+'cf_movie.sh')
            
        for pdir in proj_dirs:
            save_dir = root_dir+ tdir + mid_dirs + bdir + pdir
            if args.make_cf_histogram == 'True':
                if os.path.exists(save_dir) == False:
                    os.makedirs(save_dir)
                os.chdir(save_dir)
                if args.remove_existing_pickle != 'False':
                    subprocess.call('rm *.pkl', shell=True)
                job_string = bdir[0] + pdir.split("_")[0] + job_tag[job_it]
                nodes = n_nodes[job_it]
                if args.short_queues_only != 'False':
                    use_short = True
                else:
                    use_short = bool(np.random.randint(0,2))
                if args.astro02_only != 'False':
                    add_2 = True
                else:
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
                if queue == 'astro2_short':
                    time_str = '12:00:00'
                elif queue == 'astro_short':
                    time_str = '1-00:00:00'
                elif queue == 'astro2_long':
                    time_str = '5-00:00:00'
                else:
                    time_str = '10-00:00:00'

                if bdir.split('_')[0] == 'Bound':
                    L_bool = "True"
                else:
                    L_bool = "False"
                    
                run_string = "srun python /groups/astro/rlk/Scripts/Ramses_analysis/CF_analysis/multiplicity_fraction_res_per_bin_13.py /lustre/astro/troels/IMF_512/binary_analysis/data/ ./ -pickle cf_hist -verbose True -global_data " + global_pickles[job_it] + " -bound True -entire_sim True "
                if L_bool == "False":
                    run_string = run_string + "-lower_L 0 -upper_L 10000 -acc_lim 0 -replace_ind_type M "
                
                run_string = run_string + "-proj True -ax " + pdir.split('_')[0].lower() + " 1>cf.out00 2>&1\n"
                
                f = open('cf_pickles.sh', 'w')
                f.write("#!/bin/bash\n")
                f.write("# SLURM resource specifications\n")
                f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
                f.write("#SBATCH --job-name="+job_string+"   # shows up in the output of 'squeue'\n")
                f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
                f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
                f.write("#SBATCH --nodes="+str(nodes)+"              # number of nodes allocated for this job\n")
                f.write("#SBATCH --ntasks-per-node="+str(t_per_n)+"    # number of MPI ranks per node\n")
                f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
                f.write("#SBATCH --mail-type=ALL\n")
                
                f.write(run_string)
                f.close()

                if os.path.exists(save_dir+"cf_hist.pkl") == False:
                    subprocess.call('sbatch cf_pickles.sh', shell=True)
                    print("submitted job", save_dir+'cf_pickles.sh')
            
            if args.make_movie_frames == 'True' and os.path.exists(save_dir+'cf_pickles.sh'):
                if os.path.exists(save_dir) == False:
                    os.makedirs(save_dir)
                os.chdir(save_dir)
                if args.remove_existing_frames != 'False':
                    subprocess.call('rm movie_frame*.jpg', shell=True)
                job_string = "M" + bdir[0] + pdir.split("_")[0] + job_tag[job_it]
                nodes = n_nodes[-1*job_it]
                if args.short_queues_only != 'False':
                    use_short = True
                else:
                    use_short = bool(np.random.randint(0,2))
                if args.astro02_only != 'False':
                    add_2 = True
                else:
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
                if queue == 'astro2_short':
                    time_str = '12:00:00'
                elif queue == 'astro_short':
                    time_str = '1-00:00:00'
                elif queue == 'astro2_long':
                    time_str = '5-00:00:00'
                else:
                    time_str = '10-00:00:00'

                run_string = "srun python ~/Scripts/Ramses_analysis/CF_analysis/cf_movie.py ./ ./ -pickle cf_hist.pkl 1>cf_movie.out00 2>&1\n"
            
                f = open('cf_movie.sh', 'w')
                f.write("#!/bin/bash\n")
                f.write("# SLURM resource specifications\n")
                f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
                f.write("#SBATCH --job-name="+job_string+"   # shows up in the output of 'squeue'\n")
                f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
                f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
                f.write("#SBATCH --nodes="+str(nodes*2)+"              # number of nodes allocated for this job\n")
                f.write("#SBATCH --ntasks-per-node="+str(t_per_n)+"    # number of MPI ranks per node\n")
                f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
                f.write("#SBATCH --mail-type=ALL\n")
                
                f.write(run_string)
                f.close()

                #if os.path.exists(save_dir+"cf_hist.pkl") ==False:
                subprocess.call('sbatch cf_movie.sh', shell=True)
                print("submitted job", save_dir+'cf_movie.sh')
