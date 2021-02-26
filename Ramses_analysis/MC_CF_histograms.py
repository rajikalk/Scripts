import numpy as np
import os
import subprocess

N_samples = 20
#define priors
phi = np.arccos(np.random.random(N_samples))
theta = np.random.random(N_samples)*2*np.pi

#Calculate projection vectors
z = np.cos(phi)
xy_mag = np.sqrt(1**2 - z**2)
x = xy_mag*np.cos(theta)
y = xy_mag*np.sin(theta)
proj_vec = np.array([x, y, z]).T

for vec in proj_vec:
    save_dir = "Proj_" + str(vec[0]) + "_" + str(vec[1]) + "_" + str(vec[2])
    if os.path.exists(save_dir) == False:
        os.makedirs(save_dir)
    
    use_short = bool(np.random.randint(0,2))
    add_2 = bool(np.random.randint(0,2))
    queue = 'astro'
    if add_2:
        queue = queue + '2'
    if use_short:
        queue = queue + '_short'
    else:
        queue = queue + '_long'
    if queue != 'astro2_short':
        time_str = '1-00:00:00'
    else:
        time_str = '12:00:00'
        
    f = open(save_dir+'/cf_pickles.sh', 'w')
    f.write("#!/bin/bash\n")
    f.write("# SLURM resource specifications\n")
    f.write("# (use an extra '#' in front of SBATCH to comment-out any unused options)\n")
    f.write("#SBATCH --job-name=CFU   # shows up in the output of 'squeue'\n")
    f.write("#SBATCH --time="+time_str+"       # specify the requested wall-time\n")
    f.write("#SBATCH --partition="+queue+"  # specify the partition to run on\n")
    f.write("#SBATCH --nodes=1              # number of nodes allocated for this job\n")
    f.write("#SBATCH --ntasks-per-node=20    # number of MPI ranks per node\n")
    f.write("#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n")
    f.write("##SBATCH --exclude=node751  # avoid nodes (e.g. --exclude=node786)\n")
    f.write("#SBATCH --mail-type=ALL\n")
    proj_str='['+str(vec[0])+','+str(vec[1])+','+str(vec[2])+']'
    f.write("srun python /groups/astro/rlk/Scripts/Ramses_analysis/multiplicity_fraction_res_per_bin.py /lustre/astro/troels/IMF_512/binary_analysis/data/ "+save_dir+"/ -bound False -pickle cf_hist.pkl -proj_vec "+ proj_str +" 1>"+save_dir+"/cf.out00 2>&1\n")
    f.close()

    subprocess.call('sbatch '+save_dir+'/cf_pickles.sh', shell=True)
    print("submitted job", save_dir+'/cf_pickles.sh')
