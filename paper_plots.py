import yt
import my_module as mym
import matplotlib.pyplot as plt
#import my_fields as myf
from subprocess import call
import sys
import os
import glob
import numpy as np

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-mp", "--movie_plot", help="did you want to plot a slice of the movie?", default=True)
    parser.add_argument("-sp", "--slice_plot", help="Did you want to plot create a slice plot?", default=True)
    parser.add_argument("-pp", "--profile_plot", help="Did you want to plot a profile plot?", default=True)
    parser.add_argument("-z", "--zoom", help="Will movie be zoomed in?", default=False, type=bool)
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default=4, type=float)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", default=True)
    
    #movie plot args
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", default=True)
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-15)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-13)
    
    #slice plot args
    parser.add_argument("-sf", "--slice_field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuign volume", type=float, default=1000.)
    parser.add_argument("-dt", "--disk_thickness", help="How far above and below the midplane do you want your profile to go?", type=float, default=100.)
    parser.add_argument("-xf", "--x_field", help="x axis of the profile plot?", type=str, default="Distance_from_Center")
    parser.add_argument("-yf", "--y_field", help="y axis of the profile plot?", type=str, default="Relative_Keplerian_Velocity")
    parser.add_argument("-wf", "--weight_field", help="any weight field?", type=str, default="cell_mass")
    parser.add_argument("-log", "--logscale", help="Want to use a log scale?", type=bool, default=False)
    parser.add_argument("-pb", "--profile_bins", help="how many bins do you want for the profile?", type=int, default=100.)
    parser.add_argument("-zb", "--z_bins", help="how many z bins do you want when sampling points?", type=int, default=2.)
    parser.add_argument("-nsp", "-no_sampled_points", help="how many random points do you want to randomly sample?", type=int, default=2000)
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#======================================================================================================

args = parse_inputs()

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

if args.movie_plot == 'False':
    args.movie_plot = False
if args.movie_plot:
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/movie_script_mod.py', path, save_dir, '-pt', str(args.plot_time), '-t', args.title, '-z', str(args.zoom), '-zt', str(args.zoom_times), '-cmin', str(args.colourbar_min), '-cmax', str(args.colourbar_max), '-at', str(args.annotate_time), '-ax', 'xy', '-wr', str(args.working_rank), '-pvl', str(args.plot_velocity_legend)])

files = sorted(glob.glob(path + '*_plt_cnt*'))
file = mym.find_files([args.plot_time], files)[0]
del files
part_file = file[:-12] + 'part' + file[-5:]
ds = yt.load(file, particle_filename=part_file)
dd = ds.all_data()
movie_files = sorted(glob.glob(path + 'WIND_proj_*'))
movie_file = mym.find_files([args.plot_time], movie_files)[0]
del movie_files

if args.slice_plot == 'False':
    args.slice_plot = False
if args.slice_plot:
    save_image_name = save_dir + "Slice_Plot_time_" + str(args.plot_time) + ".eps"
    X, Y, X_vel, Y_vel = mym.initialise_grid(movie_file, args.zoom_times)
    f = h5py.File(movie_file, 'r')
    velx, vely = mym.get_quiver_arrays(f['velx_slice_xy'][:,:,0], f['vely_slice_xy'][:,:,0])
    f.close()
    ax = mym.sliceplot(ds, X, Y, args.slice_field)
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely)
    part_info = mym.get_particle_data(movie_file, axis=args.axis)
    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], 2*np.max(X), particle_mass=part_info['particle_mass'])
    plt.savefig(save_image_name, bbox_inches='tight')# , pad_inches = 0.02)
    print "created slice plot:", save_image_name


if args.profile_plot == 'False':
    args.profile_plot = False
if args.profile_plot:
    save_image_name = save_dir + "Profile_Plot_time_" + str(args.plot_time) + ".png"
    measuring_volume = ds.disk(dd['Center'], [0.0, 0.0, 1.0], (args.r_max, 'au'), (args.disk_thickness, 'au'))
    prof_x, prof_y = mym.profile_plot(measuring_volume, args.x_field, args.y_field, weight_field=args.weight_field, log=args.logscale, n_bins=args.profile_bins)
    sampled_points = mym.sampled_points(measuring_volume, args.x_field, args.y_field, bin_no=args.z_bins, no_of_points=args.no_sampled_points)
    plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
    
    plt.clf()
    fig, ax = plt.subplots()
    cm = plt.cm.get_cmap('RdYlBu')
    plot = ax.scatter(plot_array[1], plot_array[2], c=plot_array[0], alpha=0.4, cmap=cm)
    ax.plot(prof_x, prof_y[args.y_field], 'k-', linewidth=2.)
    plt.axes().set_aspect(1./plt.axes().get_data_ratio())
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('|z position| (AU)', rotation=270, labelpad=13, size=14)
    ax.set_xlabel('Cyclindral Radius (AU)', labelpad=-1)
    ax.set_ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', labelpad=-3)
    plt.axhline(y=1.0, color='k', linestyle='--')
    plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)




