import numpy as np
import pickle
import glob
from mpi4py.MPI import COMM_WORLD as CW
import yt
import os
import sys

#Set units
length_unit=yt.YTQuantity(4, 'pc')

rank = CW.Get_rank()
size = CW.Get_size()

#Generate underlying distribution:
iterations = 100
if size > 1:
    init_seed = 51092
    np.random.seed(init_seed)
    seeds = np.random.uniform(low=0.0, high=1.0, size=iterations)*10000
else:
    seeds = np.random.uniform(low=0.0, high=1.0, size=iterations)*10000
existing_dists = glob.glob('uniform_dist_*.pkl')
starting_it = len(existing_dists)
iter_range = range(starting_it, iterations)
n_stars = int(sys.argv[1])
n_bin_bound = 11
exp_bins = np.linspace(1, 6, n_bin_bound)
sep_bins = 10**exp_bins
rit = -1
ranks_per_iteration = int(sys.argv[2])
base_rank = -1

for iteration in iter_range:
    rand_pos = np.random.uniform(low=0.0, high=1.0, size=(n_stars,3))
    np.random.seed(int(seeds[iteration]))
    rand_pos = np.random.uniform(low=0.0, high=1.0, size=(n_stars,3))
    base_rank = base_rank + 1
    if (base_rank+1)*ranks_per_iteration > size:
        base_rank = 0
    if rank >= base_rank*ranks_per_iteration and rank < (base_rank+1)*ranks_per_iteration:
        #print("iteration", iteration, "going to rank", rank)
        base_rit = -1
        #separations = np.array([])
        sep_his_uni_all = np.zeros(10)
        for rand_pos_it in range(n_stars):
            base_rit = base_rit + 1
            if base_rit == ranks_per_iteration:
                base_rit = 0
            if rank == (base_rank*ranks_per_iteration + base_rit):
                diff_pos = rand_pos[rand_pos_it+1:] - rand_pos[rand_pos_it]
                update_inds = np.where(abs(diff_pos)>0.5)
                diff_pos[update_inds] = abs(diff_pos[update_inds]) - 1.0
                del update_inds
                sep = np.sqrt(diff_pos.T[0]**2 + diff_pos.T[1]**2 + diff_pos.T[2]**2)
                del diff_pos
                #separations = np.concatenate((separations,sep))
                if np.remainder(rand_pos_it,100) == 0:
                    print("updated separation for star pos", rand_pos_it, "for iteration", iteration, "on rank", rank)
                sep_his_uni, bins = np.histogram(sep*length_unit.in_units('AU'), bins=sep_bins)
                del sep
                del bins
                sep_his_uni_all = sep_his_uni_all + sep_his_uni
                del sep_his_uni
        
        pickle_name = 'seps_iter_'+str(iteration)+'_rank_'+("%06d" % rank)+'.pkl'
        file = open(pickle_name, 'wb')
        pickle.dump((sep_his_uni_all), file)
        file.close()
        del sep_his_uni_all
        print("Saved separations for No", iteration, "of", iterations, "on rank", rank)
        
        sys.stdout.flush()
        CW.Barrier()
        
        if rank == base_rank*ranks_per_iteration:
            pickle_files = sorted(glob.glob('seps_iter_'+str(iteration)+'_rank_*.pkl'))
            sep_his_uni_all = np.zeros(10)
            for pickle_file in pickle_files:
                file_open = open(pickle_file, 'rb')
                sep_his_uni = pickle.load(file_open)
                file_open.close()
                sep_his_uni_all = sep_his_uni_all + sep_his_uni
                del sep_his_uni
                os.remove(pickle_file)
                print('read', pickle_file)
            
            pickle_name = 'uniform_dist_'+("%06d" % iteration)+'.pkl'
            file = open(pickle_name, 'wb')
            pickle.dump((sep_his_uni_all), file)
            file.close()
            del pickle_name
            del sep_his_uni_all
            #RR = sep_his_uni/(n_stars*(n_stars-1))
            print("created histogram No", iteration, "of", iterations)

sys.stdout.flush()
CW.Barrier()

if rank == 0:
    uni_hists = []
    pickle_files = glob.glob('uniform_dist_*.pkl')
    for pickle_file in pickle_files:
        file_open = open(pickle_file, 'rb')
        sep_his_uni = pickle.load(file_open)
        file_open.close()
        uni_hists.append(sep_his_uni)
        os.remove(pickle_file)
    
    pickle_name = 'all_uniform_dist.pkl'
    file = open(pickle_name, 'wb')
    pickle.dump((uni_hists), file)
    file.close()

sys.stdout.flush()
CW.Barrier()
print('Finished deriving underlying distribution')
