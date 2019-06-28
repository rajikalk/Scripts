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
    parser.add_argument("-c", "--center", help="What is the center?", default=0, type=int)
    parser.add_argument("-a", "--semimajor_axis", help="what is the semimajor axis of the binary?", default=0, type=float)
    parser.add_argument("-mr", "--max_r", help="what is the max radius?", default=None, type=float)
    parser.add_argument("-bins", "--no_of_bins", help="number of bins to use", default=None, type=int)
    parser.add_argument("-ab", "--adaptive_bins", help="do you want to use adaptive bin sizes?", default=True, type=str)
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    args = parse_inputs()
    center = int(args.center)
    if args.adaptive_bins == 'False':
        args.adaptive_bins = False
    a = args.semimajor_axis
    max_radius = args.max_r
    
    file = open(args.file, 'r')
    loaded_fields = pickle.load(file)
    distance = loaded_fields[0]
    cell_mass = loaded_fields[1]
    part_mass = loaded_fields[2]
    
    rank = CW.Get_rank()
    size = CW.Get_size()

    dist_min = np.min(distance)
    dist_max = np.max(distance)
    if args.adaptive_bins:
        rs = [0]
        rs = rs + list(set(distance[distance<=max_radius]))
        rs = np.sort(rs)
        bin_freq = len(rs)/500
        if bin_freq > 0:
            rs = rs[::bin_freq]
        rs = np.array(rs)
        rs = np.append(rs, max_radius)
    else:
        rs = np.linspace(0.0, max_radius, args.no_of_bins)
    bin_size = rs[-1] - rs[-2]
    #print "bin_size =", bin_size/1.49597870751e+13
    #print "max_radius =", max_radius/1.49597870751e+13
    rs = np.append(rs, rs[-1]+bin_size)
    enclosed_mass = np.array(np.zeros(np.shape(distance)))

    rit = 1
    printed = False
    print_cen = False
    for r in range(len(rs))[1:]:
        if rank == rit:
            enclosed_mass = np.array(np.zeros(np.shape(distance)))
            ind = np.where((distance >= rs[r-1]) & (distance < rs[r]))[0]
            enclosed_dist = np.where(distance < rs[r-1])[0]
            if len(enclosed_dist) == 0:
                enclosed_dist = np.where(distance < (dist_min+1))
            enclosed_mass_val = np.sum(cell_mass[enclosed_dist])
            if center != 0:
                enclosed_mass_val = enclosed_mass_val + part_mass[center-1]
                if print_cen == False:
                    #print "Centered on particle with mass", part_mass[center-1]/1.98841586e+33
                    print_cen = True
            if center !=0 and rs[r] > a and len(part_mass)>1:
                if center == 1:
                    enclosed_mass_val = enclosed_mass_val + part_mass[1]
                else:
                    enclosed_mass_val = enclosed_mass_val + part_mass[0]
                if printed == False:
                    #print "Added other particle with mass", part_mass[0]/1.98841586e+33
                    printed = True
            elif center == 0 and rs[r] > a/2. and len(part_mass)>0:
                enclosed_mass_val = enclosed_mass_val + np.sum(part_mass)
                if printed == False:
                    #print "Added both particles with mass", np.sum(part_mass)/1.98841586e+33
                    printed = True
            enclosed_mass[ind] = enclosed_mass_val
            #print "enclosed mass =", enclosed_mass_val/1.98841586e+33, ", Radius =", rs[r]/1.49597870751e+13, "on rank", rank
            CW.send(enclosed_mass, dest=0, tag=rank)
        if rank == 0:
            enclosed_mass_add = CW.recv(source=rit, tag=rit)
            enclosed_mass = enclosed_mass + enclosed_mass_add
        rit = rit + 1
        if rit == size:
            rit = 1

    if rank == 0:
        os.remove(args.file)
        file = open(args.save_file, 'w+')
        print("pickle file:", args.save_file)
        pickle.dump(enclosed_mass, file)
        file.close()

if __name__ == '__main__': main()
