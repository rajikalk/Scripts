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
    parser.add_argument("-bs", "--no_of_bins", help="how many bins", default=None, type=float)
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    args = parse_inputs()
    no_of_bins = args.no_of_bins
    
    file = open(args.file, 'r')
    z, B= pickle.load(file)
    
    dz = list(set(z))
    rs = [0]
    rs = rs + list(np.sort(dz))

    rank = CW.Get_rank()
    size = CW.Get_size()

    max_height = np.max(z)
    #rs = np.linspace(np.min(z), max_height, no_of_bins)
    bin_size = rs[-1] - rs[-2]
    print "bin_size =", bin_size/1.49597870751e+13
    print "max_radius =", max_height/1.49597870751e+13
    rs = rs.append(rs[-1]+bin_size)
    rs = np.array(rs)
    gradient = np.array(np.zeros(np.shape(z)))

    rit = 1
    for r in range(len(rs))[1:]:
        if rank == rit:
            grad = np.array(np.zeros(np.shape(z)))
            if r == 1:
                B_0 = 0.0
            else:
                shell_0 = np.where((z >= rs[r-2])&(z < rs[r-1]))[0]
                if len(shell_0) == 0:
                    B_0 = 0.0
                else:
                    B_0 = np.mean(B[shell_0])
            shell_1 = np.where((z >= rs[r-1])&(z < rs[r]))[0]
            if len(shell_1) == 0:
                B_1 = 0.0
            else:
                B_1 = np.mean(B[shell_1])
            if rs[r] == rs[-1]:
                B_2 = 0.0
            else:
                shell_2 = np.where((z >= rs[r])&(z < rs[r+1]))[0]
                if len(shell_2) == 0:
                    B_2 = 0.0
                else:
                    B_2 = np.mean(B[shell_2])
            print "B_0, B_1, B_2 =", B_0, B_1, B_2
            if rs[r] > 0.0:
                B_grad = (((B_1 - B_0)/(rs[r] - rs[r-1])) + ((B_2 - B_1)/(rs[r+1] - rs[r])))/2.
            else:
                B_grad = (((B_0 - B_1)/(rs[r-1]-rs[r])) + ((B_1 - B_2)/(rs[r]-rs[r+1])))/2.
            print "B_grad =", B_grad, "on rank", rank
            grad[shell_1] = B_grad
            CW.send(grad, dest=0, tag=rank)
        if rank == 0:
            gradient_add = CW.recv(source=rit, tag=rit)
            gradient = gradient + gradient_add
        rit = rit + 1
        if rit == size:
            rit = 1

    if rank == 0:
        os.remove(args.file)
        file = open(args.save_file, 'w+')
        pickle.dump(gradient, file)
        file.close()

if __name__ == '__main__': main()
