#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
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
    parser.add_argument("-thickness", "--slice_thickness", help="How thick would you like your yt_projections to be? default 100AU", type=float, default=100.)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
def main():
    
    file = sys.argv[1]
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    L = [0.0, 0.0, 1.0]
    myf.set_normal(L)
    center_pos = np.array([0.0, 0.0, 0.0])

    args = parse_inputs()
    field = args.field
    dd = ds.all_data()
    center_vel = dd['Center_Velocity']
    center_pos = dd['Center_Position']
    part_pos = dd['All_Particle_Positions']
    part_mass = dd['All_Particle_Masses']
    print "center vel =", myf.get_center_vel()
    print "center pos for field calculations", myf.get_center_pos()
    print "particle positions are:", myf.get_part_pos()
    print "particle masses are:", myf.get_part_mass()
    center_pos = np.array([0.0, 0.0, 0.0])
    proj = yt.OffAxisProjectionPlot(ds, L, [field, 'cell_mass', 'velx_mw', 'vely_mw', 'magx_mw', 'magy_mw'], center=(center_pos, 'AU'), width=(args.ax_lim, 'AU'), depth=(args.slice_thickness, 'AU'))
    image = (proj.frb.data[simfo['field']]/thickness.in_units('cm')).value
    velx_full = (proj.frb.data[('gas', 'Projected_Velocity_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
    vely_full = (proj.frb.data[('gas', 'velz_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).value
    magx = (proj.frb.data[('gas', 'Projected_Magnetic_Field_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
    magy = (proj.frb.data[('gas', 'magz_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).value
    mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm')).value
    proj.save()



if __name__ == '__main__': main()
