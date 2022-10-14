#!/usr/bin/env python
import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="z")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
args = parse_inputs()

#Get movie files
movie_files = sorted(glob.glob(input_dir + '*_plt_cnt*'))

#Now let's iterate over the files and get the images we want to plot
for fn in yt.parallel_objects(movie_files, njobs=int(size/5)):
    part_file = 'part'.join(fn.split('plt_cnt'))
    ds = yt.load(fn, particle_filename=part_file)
    
    #make list of projection fields: density, velocity, magnetic field
    proj_field_list = [('flash', 'dens')] + \
        [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+args.axis not in field[1])] + \
        [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+args.axis not in field[1])]
    
    #This is the dictionary where the projected arrays will be saved:
    proj_dict = {}
    for field in proj_field_list:
        proj_dict.update({field[1]:[]})
    
    for field in yt.parallel_objects(proj_field_list):
        proj = yt.ProjectionPlot(ds, args.axis, , method='integrate')
        thickness = (proj.bounds[1] - proj.bounds[0]).in_cgs() #MIGHT HAVE TO UPDATE THIS LATER
        proj_array = proj.frb.data[field].in_cgs()/thickness
        import pdb
        pdb.set_trace()
