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
    parser.add_argument("-bs", "--no_of_bins", help="how many bins", default=100, type=float)
    parser.add_argument("-co", "--coordinates", help="What is the coordinates?", default='cylindrical', type=str)
    parser.add_argument("-a", "--semimajor_axis", help="what is the semimajor axis of the binary?", default=0, type=float)
    parser.add_argument("-mr", "--max_r", help="what is the max radius?", default=None, type=float)
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    args = parse_inputs()
    center = int(args.center)
    no_of_bins = args.no_of_bins
    coordinates = args.coordinates
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
    rs = np.linspace(0.0, max_radius, no_of_bins)
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
                    print "Centered on particle with mass", part_mass[center-1]/1.98841586e+33
                    print_cen = True
            if center !=0 and rs[r] > a:
                if center == 1:
                    enclosed_mass_val = enclosed_mass_val + part_mass[1]
                else:
                    enclosed_mass_val = enclosed_mass_val + part_mass[0]
                if printed == False:
                    print "Added other particle with mass", part_mass[0]/1.98841586e+33
                    printed = True
            elif center == 0 and rs[r] > a/2. and len(part_mass)>0:
                enclosed_mass_val = enclosed_mass_val + np.sum(part_mass)
                if printed == False:
                    print "Added both particles with mass", np.sum(part_mass)/1.98841586e+33
                    printed = True
            enclosed_mass[ind] = enclosed_mass_val
            print "enclosed mass =", enclosed_mass_val/1.98841586e+33, ", Radius =", rs[r]/1.49597870751e+13, "on rank", rank
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
        print "pickle file:", args.save_file
        pickle.dump(enclosed_mass, file)
        file.close()

if __name__ == '__main__': main()
