#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import argparse
import pickle
import os

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="file to use")
    parser.add_argument("-sf", "--save_file", help="file to save to")
    parser.add_argument("-mr", "--max_r", help="what is the max radius?", default=None, type=float)
    parser.add_argument("-bins", "--no_of_bins", help="number of bins to use", default=None, type=int)
    parser.add_argument("-ab", "--adaptive_bins", help="do you want to use adaptive bin sizes?", default=True, type=str)
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    args = parse_inputs()
    if args.adaptive_bins == 'False':
        args.adaptive_bins = False
    max_radius = args.max_r
    
    file = open(args.file, 'r')
    loaded_fields = pickle.load(file)
    distance = loaded_fields[0]
    y = loaded_fields[1]
    bin_data = loaded_fields[2]
    
    rank = CW.Get_rank()
    size = CW.Get_size()

    dist_min = np.min(distance)
    dist_max = np.max(distance)
    if args.adaptive_bins:
        rs = [0] + list(set(distance[bin_data<=max_radius]))
        rs = np.sort(rs)
        rs = rs[::2]
        rs = np.array(rs)
        rs = np.append(rs, max_radius)
    else:
        rs = np.linspace(0.0, max_radius, args.no_of_bins)
    bin_size = rs[-1] - rs[-2]
    rs = np.append(rs, rs[-1]+bin_size)
    gradient = np.array(np.zeros(np.shape(distance)))
    #print "RS:", rs

    rit = 1
    printed = False
    print_cen = False
    for r in range(len(rs)):
        if rank == rit:
            grad_add = np.array(np.zeros(np.shape(distance)))
            if r-1 < 0:
                r_0 = rs[r]
            else:
                r_0 = rs[r-1]
            r_1 = rs[r]
            if r+1 == len(rs):
                r_2 = rs[r]
            else:
                r_2 = rs[r+1]
            if r+2 == len(rs)+1:
                r_3 = rs[r]
            elif r+2 == len(rs):
                r_3 = rs[r+1]
            else:
                r_3 = rs[r+2]
            mid_01 = (r_1 + r_0)/2.
            mid_23 = (r_3 + r_2)/2.
            shell_01 = np.where((distance >= r_0) & (distance < r_1))[0]
            shell_12 = np.where((distance >= r_1) & (distance < r_2))[0]
            shell_23 = np.where((distance >= r_2) & (distance < r_3))[0]
            if len(shell_01) == 0:
                y_01 = 0.0
            else:
                y_01 = np.mean(y[shell_01])
            if len(shell_23) == 0:
                y_23= 0.0
            else:
                y_23 = np.mean(y[shell_23])
            grad_val = (y_23 - y_01)/(2.*(mid_23 - mid_01))
            #if rank == 1:
            #print "r_0, r_1, r_2, r_3:", r_0, r_1, r_2, r_3
            #print "mid_01, mid_12, mid_23:", mid_01, mid_12, mid_23
            #print "y_01, y_12, y_23:", y_01, y_12, y_23, "on rank", rank
            #print "grad_1, grad_2, average", grad_1, grad_2, grad_val, "on rank", rank
            #print "Gradient =", grad_val, "at Distance =", np.mean([mid_01, mid_23]), "on rank", rank
            grad_add[shell_12] = grad_val
            #grad_add[shell] = grad_val
            CW.send(grad_add, dest=0, tag=rank)
        if rank == 0:
            grad_add = CW.recv(source=rit, tag=rit)
            gradient = gradient + grad_add
        rit = rit + 1
        if rit == size:
            rit = 1

    if rank == 0:
        os.remove(args.file)
        file = open(args.save_file, 'w+')
        print "pickle file:", args.save_file
        pickle.dump(gradient, file)
        file.close()

if __name__ == '__main__': main()
