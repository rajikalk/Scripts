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
    y = loaded_fields[1]
    
    rank = CW.Get_rank()
    size = CW.Get_size()

    dist_min = np.min(distance)
    dist_max = np.max(distance)
    if args.adaptive_bins:
        rs = [0]
        rs = rs + list(set(distance[distance<=max_radius]))
        rs = np.sort(rs)
        rs = rs[::2]
        rs = np.array(rs)
        rs = np.append(rs, max_radius)
    else:
        rs = np.linspace(0.0, max_radius, args.no_of_bins)
    bin_size = rs[-1] - rs[-2]
    rs = np.append(rs, rs[-1]+bin_size)
    gradient = np.array(np.zeros(np.shape(distance)))

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
            shell_0_1 = np.where((distance >= r_0) & (distance < r_1))[0]
            shell_1_2 = np.where((distance >= r_1) & (distance < r_2))[0]
            shell_2_3 = np.where((distance >= r_2) & (distance < r_3))[0]
            if len(shell_0_1) == 0:
                y_0_1 = 0.0
            else:
                y_0_1 = np.mean(y[shell_0_1])
            if len(shell_1_2) == 0:
                y_1_2= 0.0
            else:
                y_1_2 = np.mean(y[shell_1_2])
            if len(shell_2_3) == 0:
                y_2_3= 0.0
            else:
                y_2_3 = np.mean(y[shell_2_3])
            grad_val = ((y_1_2 - y_0_1)/(r_1 - r_0) + (y_2_3 - y_1_2)/(r_2 - r_1))/2.
            print("Gradient =", grad_val, "at Distance =", r_1, "on rank", rank)
            grad_add[shell_1_2] = grad_val
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
            print("enclosed mass =", enclosed_mass_val/1.98841586e+33, ", Radius =", rs[r]/1.49597870751e+13, "on rank", rank)
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
