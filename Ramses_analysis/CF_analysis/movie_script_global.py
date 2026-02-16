#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
#import my_ramses_fields_short as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import gc

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-22')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-19)
    parser.add_argument("-axis", "--perp_axis", help="What axis do you want to take the projection along?", type=str, default="z")
    parser.add_argument("-make_proj", "--make_projection", help="do you want to make projections or just skip onto the particles?", type=str, default="True")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)
    
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
global_pickle = sys.argv[3]
if os.path.exists(save_dir) == False and rank==0:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

#Set some plot variables independant on data files
#File files
sink_files = sorted(glob.glob(input_dir+"output*/*.dat"))
files = sorted(glob.glob(input_dir+"*/info*.txt"))
rm_files = []
for info_name in files:
    sink_file = info_name.split('info')[0]+'stars_output.dat'
    if sink_file not in sink_files:
        rm_files.append(info_name)
for rm_file in rm_files:
    files.remove(rm_file)
del sink_files
gc.collect()

print('Collected files on rank', rank)

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = global_pickle.split('/G')[-1].split('/')[0] #save_dir.split('G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    
del simulation_density_id
gc.collect()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

print('Calculated the units on rank', rank)
sys.stdout.flush()
CW.Barrier()

#load global sink data

usable_files = files
del files
gc.collect()
sys.stdout.flush()
CW.Barrier()
'''
rit = -1
for fn_it in range(len(usable_files)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
'''
if args.make_projection == "True":
    for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(100))):
        pickle_file = save_dir + "movie_frame_" + ("%06d" % fn_it) + "_proj.pkl"
        if os.path.exists(pickle_file) == False:
            fn = usable_files[fn_it]
            ds = yt.load(fn, units_override=units_override)
            
            time_real = ds.current_time.in_units('yr')
            time_val = np.round(time_real.in_units('yr'))
            
            proj = yt.ProjectionPlot(ds, args.perp_axis, ("ramses", "Density"), width=units['length_unit'], method='integrate')
            proj_array = np.array(proj.frb.data[("ramses", "Density")]/units['length_unit'].in_units('cm'))
            image = proj_array*units['density_unit'].in_units('g/cm**3')
            print('Made projection of file', fn, 'on rank', rank)
            del proj
            del proj_array
            gc.collect()
            '''
            if Grho == 200 or Grho == 400:
                del dd
                file_open = open(global_pickle, 'rb')
                try:
                    global_data = pickle.load(file_open,encoding="latin1")
                except:
                    file_open.close()
                    import pickle5 as pickle
                    file_open = open(global_pickle, 'rb')
                    global_data = pickle.load(file_open,encoding="latin1")
                file_open.close()
                time_arr = global_data['time'].T[0]*scale_t.in_units('yr')
                global_t_ind = np.argmin(abs(time_arr - time_real))
                n_stars = np.where(global_data['m'][global_t_ind]>0)[0]
                particle_positions = np.array([global_data['x'][global_t_ind][n_stars], global_data['y'][global_t_ind][n_stars], global_data['z'][global_t_ind][n_stars]])*scale_l
                particle_masses = global_data['m'][global_t_ind][n_stars]*units['mass_unit']
                del global_t_ind
                del n_stars
            
            dd = ds.all_data()
            particle_x_pos = dd['sink_particle_posx']
            particle_y_pos = dd['sink_particle_posy']
            #particle_masses = dd['sink_particle_mass']
            del ds
            del dd
            '''
            del time_real
            gc.collect()
            

            if np.remainder(rank,100) == 0:
                file = open(pickle_file, 'wb')
                #pickle.dump((image, time_val, particle_positions, particle_masses), file)
                pickle.dump((image, time_val), file)
                file.close()
                print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
            del image
            #del particle_x_pos
            #del particle_y_pos
            #del particle_masses
            del time_val
            gc.collect()
                
print('FINISHED MAKING YT PROJECTIONS ON RANK', rank)

sys.stdout.flush()
CW.Barrier()

from pyramses import rsink
'''
sys.stdout.flush()
CW.Barrier()
if args.perp_axis == 'z' or args.perp_axis == 'y':
    def _sink_particle_posx(field, data):
        """
        Retrieve particle x position from .snktxt file
        """
        particle_posx = []
        if np.shape(data['x']) == (16, 16, 16):
            particle_posx = yt.YTArray(np.array(particle_posx), "pc")
        else:
            file_no = int(data.ds.directory.split('output_')[-1])
            datadir = data.ds.directory.split('output_')[0]
            loaded_sink_data = rsink(file_no, datadir=datadir)
            del file_no
            del datadir
            gc.collect()
            """
            xlim= [np.min(data['x']-data['dx']).in_units('pc').value, np.max(data['x']+data['dx']).in_units('pc').value]
            ylim= [np.min(data['y']-data['dy']).in_units('pc').value, np.max(data['y']+data['dy']).in_units('pc').value]
            zlim= [np.min(data['z']-data['dz']).in_units('pc').value, np.max(data['z']+data['dz']).in_units('pc').value]
            posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
            posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
            posz = loaded_sink_data['z']*data.ds.length_unit.in_units("pc").value
            del loaded_sink_data
            gc.collect()
            xinds = np.argwhere((posx>xlim[0])&(posx<xlim[1]))
            yinds = np.argwhere((posy>ylim[0])&(posy<ylim[1]))
            zinds = np.argwhere((posz>zlim[0])&(posz<zlim[1]))
            del xlim
            del ylim
            del zlim
            del posy
            del posz
            gc.collect()
            usable_inds = set(xinds).intersection(set(yinds).intersection(set(zinds)))
            del posx
            gc.collect()
            particle_posx = yt.YTArray(np.array(posx[usable_inds]), "pc")
            """
            particle_posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
            particle_posx = yt.YTArray(particle_posx, "pc")
        return particle_posx

    yt.add_field("sink_particle_posx", function=_sink_particle_posx, units=r"pc")

sys.stdout.flush()
CW.Barrier()

if args.perp_axis == 'z' or args.perp_axis == 'x':
    def _sink_particle_posy(field, data):
        """
        Retrieve particle y position from .snktxt file
        """
        particle_posy = []
        if np.shape(data['y']) == (16, 16, 16):
            particle_posy = yt.YTArray(np.array(particle_posy), "pc")
        else:
            file_no = int(data.ds.directory.split('output_')[-1])
            datadir = data.ds.directory.split('output_')[0]
            loaded_sink_data = rsink(file_no, datadir=datadir)
            del file_no
            del datadir
            gc.collect()
            """
            xlim= [np.min(data['x']-data['dx']).in_units('pc').value, np.max(data['x']+data['dx']).in_units('pc').value]
            ylim= [np.min(data['y']-data['dy']).in_units('pc').value, np.max(data['y']+data['dy']).in_units('pc').value]
            zlim= [np.min(data['z']-data['dz']).in_units('pc').value, np.max(data['z']+data['dz']).in_units('pc').value]
            posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
            posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
            posz = loaded_sink_data['z']*data.ds.length_unit.in_units("pc").value
            del loaded_sink_data
            gc.collect()
            xinds = np.argwhere((posx>xlim[0])&(posx<xlim[1]))
            yinds = np.argwhere((posy>ylim[0])&(posy<ylim[1]))
            zinds = np.argwhere((posz>zlim[0])&(posz<zlim[1]))
            del xlim
            del ylim
            del zlim
            del posx
            del posz
            gc.collect()
            usable_inds = set(xinds).intersection(set(yinds).intersection(set(zinds)))
            del posy
            gc.collect()
            particle_posy = yt.YTArray(np.array(posy[usable_inds]), "pc")
            """
            particle_posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
            particle_posy = yt.YTArray(particle_posy, "pc")
        return particle_posy

    yt.add_field("sink_particle_posy", function=_sink_particle_posy, units=r"pc")

sys.stdout.flush()
CW.Barrier()

if args.perp_axis == 'y' or args.perp_axis == 'x':
    def _sink_particle_posz(field, data):
        """
        Retrieve particle y position from .snktxt file
        """
        particle_posz = []
        if np.shape(data['z']) == (16, 16, 16):
            particle_posz = yt.YTArray(np.array(particle_posz), "pc")
        else:
            file_no = int(data.ds.directory.split('output_')[-1])
            datadir = data.ds.directory.split('output_')[0]
            loaded_sink_data = rsink(file_no, datadir=datadir)
            del file_no
            del datadir
            gc.collect()
            """
            xlim= [np.min(data['x']-data['dx']).in_units('pc').value, np.max(data['x']+data['dx']).in_units('pc').value]
            ylim= [np.min(data['y']-data['dy']).in_units('pc').value, np.max(data['y']+data['dy']).in_units('pc').value]
            zlim= [np.min(data['z']-data['dz']).in_units('pc').value, np.max(data['z']+data['dz']).in_units('pc').value]
            posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
            posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
            posz = loaded_sink_data['z']*data.ds.length_unit.in_units("pc").value
            del loaded_sink_data
            gc.collect()
            xinds = np.argwhere((posx>xlim[0])&(posx<xlim[1]))
            yinds = np.argwhere((posy>ylim[0])&(posy<ylim[1]))
            zinds = np.argwhere((posz>zlim[0])&(posz<zlim[1]))
            del xlim
            del ylim
            del zlim
            del posx
            del posy
            gc.collect()
            usable_inds = set(xinds).intersection(set(yinds).intersection(set(zinds)))
            del posz
            gc.collect()
            particle_posz = yt.YTArray(np.array(posz[usable_inds]), "pc")
            """
            particle_posz = loaded_sink_data['z']*data.ds.length_unit.in_units("pc").value
            particle_posz = yt.YTArray(particle_posz, "pc")
        return particle_posz

    yt.add_field("sink_particle_posz", function=_sink_particle_posz, units=r"pc")
    
def _sink_particle_mass(field, data):
    """
    Retrieve particle mass from .snktxt file
    """
    particle_mass = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_mass = yt.YTArray(np.array(particle_mass), "Msun")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        del file_no
        del datadir
        gc.collect()
        particle_mass = loaded_sink_data['m']*data.ds.mass_unit.in_units("Msun").value
        particle_mass = yt.YTArray(np.array(particle_mass), "Msun")
    return particle_mass

yt.add_field("sink_particle_mass", function=_sink_particle_mass, units=r"Msun")
'''
sys.stdout.flush()
CW.Barrier()
gc.collect()

#for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(48))):
rit = -1
for fn_it in range(len(usable_files)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        pickle_file = save_dir + "movie_frame_" + ("%06d" % fn_it) + "_part.pkl"
        if os.path.exists(pickle_file) == False:
            fn = usable_files[fn_it]
            file_no = int(fn.split('output_')[-1].split('/')[0])
            datadir = fn.split('output_')[0]
            loaded_sink_data = rsink(file_no, datadir=datadir)
            particle_masses = loaded_sink_data['m']*units['mass_unit']
            if args.perp_axis == "z":
                particle_x_pos = loaded_sink_data['x']*units['length_unit']
                particle_y_pos = loaded_sink_data['y']*units['length_unit']
            elif args.perp_axis == "y":
                particle_x_pos = loaded_sink_data['x']*units['length_unit']
                particle_y_pos = loaded_sink_data['z']*units['length_unit']
            else:
                particle_x_pos = loaded_sink_data['y']*units['length_unit']
                particle_y_pos = loaded_sink_data['z']*units['length_unit']
            gc.collect()
            #particle_masses = dd['sink_particle_mass']

            #if np.remainder(rank,48) == 0:
            file = open(pickle_file, 'wb')
            #pickle.dump((image, time_val, particle_positions, particle_masses), file)
            pickle.dump((particle_x_pos, particle_y_pos, particle_masses), file)
            file.close()
            print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
            del particle_x_pos
            del particle_y_pos
            gc.collect()
        
print('FINISHED GETTING PARTICLE DATA ON RANK', rank)

sys.stdout.flush()
CW.Barrier()

proj_pickles = sorted(glob.glob(save_dir + "movie_frame_*_proj.pkl"))
part_pickles = sorted(glob.glob(save_dir + "movie_frame_*_part.pkl"))

x = np.linspace(-1*units['length_unit']/2., units['length_unit']/2., 800)
y = np.linspace(-1*units['length_unit']/2., units['length_unit']/2., 800)
xlim = [x[0], x[-1]]
ylim = [y[0], y[-1]]
X, Y = np.meshgrid(x, y)

cbar_max = args.colourbar_max
try:
    cbar_min = float(args.colourbar_min)
except:
    cbar_min = float(args.colourbar_min[1:])

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import transforms
import matplotlib.patheffects as path_effects
from matplotlib import ticker

rit = -1
for pick_it in range(len(proj_pickles)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        save_name = save_dir + "movie_frame_" + ("%06d" % pick_it)
        if os.path.exists(save_name+".jpg") == False:
            file = open(proj_pickles[pick_it], 'rb')
            image, time_val = pickle.load(file)
            file.close()
            
            plt.clf()
            fig, ax = plt.subplots()
            if args.perp_axis == "z":
                ax.set_xlabel("$x$ (AU)", labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel("$y$ (AU)", fontsize=args.text_font) #, labelpad=-20
            elif args.perp_axis == "y":
                ax.set_xlabel("$x$ (AU)", labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel("$z$ (AU)", fontsize=args.text_font) #, labelpad=-20
                image = image.T
            else:
                ax.set_xlabel("$y$ (AU)", labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel("$z$ (AU)", fontsize=args.text_font) #, labelpad=-20
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar(plot, pad=0.0)
            
            file = open(part_pickles[pick_it], 'rb')
            particle_x_pos, particle_y_pos, particle_masses = pickle.load(file)
            file.close()
            
            SFE_val = np.sum(particle_masses)/units['mass_unit']
            SFE_percent = np.round(SFE_val*100, decimals=1)
            if SFE_percent<=5:
                ax.scatter((particle_x_pos.value - 2), (particle_y_pos.value - 2), color='c', s=0.5)
                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)

                time_string = "$SFE$="+str(SFE_percent)+"\%"
                time_string_raw = r"{}".format(time_string)
                time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
                
                file = open(proj_pickles[pick_it].split('_proj')[0]+'.pkl', 'wb')
                pickle.dump((X, Y, image, particle_x_pos, particle_y_pos), file)
                file.close()
                
                plt.savefig(save_name + ".jpg", format='jpg', bbox_inches='tight')
                plt.savefig(save_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame on rank', rank, 'at time of', str(time_val), 'to save_dir:', save_name + '.jpg')
            
sys.stdout.flush()
CW.Barrier()
print("Finished making frames")
