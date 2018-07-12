from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import pickle
import sys
import yt
yt.enable_parallelism()
import my_fields as myf

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="")
    parser.add_argument("-o", "--output_pickle_name", help="What will you save your output files as?", default="projection.pkl")
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number? (in AU)", type=float, default=300.)
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 100AU", type=float, default=300.)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    rank = CW.Get_rank()
    size = CW.Get_size()
    #comm = MPI.COMM_WORLD
    
    # Read in directories:
    file = sys.argv[1]
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    L = [0.0, 0.0, 1.0]
    center_pos = np.array([0.0, 0.0, 0.0])

    args = parse_inputs()
    field = args.field
    x_width = args.ax_lim
    thickness = args.slice_thickness
    print "reading in fields"
    fields = ds.derived_field_list
    print "Loaded all data"
    dd = ds.all_data()
    print "calculating field"
    proj_field = dd[field]
    print "creating projection of field", field
    proj = yt.OffAxisProjectionPlot(ds, L, [field, 'cell_mass', 'velx_mw', 'vely_mw', 'magx_mw', 'magy_mw'], center=(center_pos, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))

    image = (proj.frb.data[field]/thickness.in_units('cm')).value
    velx_full = (proj.frb.data[('gas', 'velx_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
    vely_full = (proj.frb.data[('gas', 'vely_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
    magx = (proj.frb.data[('gas', 'magx_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
    magy = (proj.frb.data[('gas', 'magy_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
    mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm')).value
    velx_full = velx_full/mass
    vely_full = vely_full/mass
    magx = magx/mass
    magy = magy/mass
    del mass
    
    if rank == 0:
        pickle_file = args.output_pickle_name
        file = open(pickle_file, 'w+')
        pickle.dump((image, velx_full, vely_full, magx, magy), file)
        file.close()
        print "Created Projection Pickle:", pickle_file, "for  file:", file

if __name__ == '__main__': main()
