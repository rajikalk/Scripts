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



#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

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

simulation_density_id = save_dir.split('G')[-1].split('/')[0]

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
for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(100))):
    pickle_file = save_dir + "movie_frame_" + ("%06d" % fn_it) + "_proj.pkl"
    if os.path.exists(pickle_file) == False:
        fn = usable_files[fn_it]
        ds = yt.load(fn, units_override=units_override)
        
        time_real = ds.current_time.in_units('yr')
        time_val = np.round(time_real.in_units('yr'))
        
        proj = yt.ProjectionPlot(ds, "z", ("ramses", "Density"), width=units['length_unit'], method='integrate')
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

sys.stdout.flush()
CW.Barrier()

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
        particle_posx = loaded_sink_data['x']*data.ds.length_unit.in_units("pc").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_posx.append(float(row[1])*data.ds.length_unit.value)
                except:
                    continue
        '''
        particle_posx = yt.YTArray(np.array(particle_posx), "pc")
    return particle_posx

yt.add_field("sink_particle_posx", function=_sink_particle_posx, units=r"pc")

sys.stdout.flush()
CW.Barrier()

def _sink_particle_posy(field, data):
    """
    Retrieve particle y position from .snktxt file
    """
    particle_posy = []
    if np.shape(data['x']) == (16, 16, 16):
        particle_posy = yt.YTArray(np.array(particle_posy), "pc")
    else:
        file_no = int(data.ds.directory.split('output_')[-1])
        datadir = data.ds.directory.split('output_')[0]
        loaded_sink_data = rsink(file_no, datadir=datadir)
        del file_no
        del datadir
        gc.collect()
        particle_posy = loaded_sink_data['y']*data.ds.length_unit.in_units("pc").value
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_posy.append(float(row[2])*data.ds.length_unit.value)
                except:
                    continue
        '''
        particle_posy = yt.YTArray(np.array(particle_posy), "pc")
    return particle_posy

yt.add_field("sink_particle_posy", function=_sink_particle_posy, units=r"pc")

sys.stdout.flush()
CW.Barrier()

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
        '''
        particle_file = glob.glob(data.ds.directory + '/*.snktxt')
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
        time_val = None
        with open(particle_file[0], 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if time_val == None:
                    time_val = yt.YTQuantity(float(row[1])*data.ds.time_unit.value, 's')
                    continue
                try:
                    dummy = float(row[0])
                    particle_mass.append(float(row[7])*data.ds.mass_unit.value)
                except:
                    continue
        '''
        particle_mass = yt.YTArray(np.array(particle_mass), "Msun")
    return particle_mass

yt.add_field("sink_particle_mass", function=_sink_particle_mass, units=r"Msun")

sys.stdout.flush()
CW.Barrier()

#for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/(48))):
for fn_it in range(len(usable_files)):
    pickle_file = save_dir + "movie_frame_" + ("%06d" % fn_it) + "_part.pkl"
    if os.path.exists(pickle_file) == False:
        fn = usable_files[fn_it]
        ds = yt.load(fn, units_override=units_override)
        dd = ds.all_data()
        particle_x_pos = dd['sink_particle_posx']
        particle_y_pos = dd['sink_particle_posy']
        del ds
        del dd
        gc.collect()
        #particle_masses = dd['sink_particle_mass']

        #if np.remainder(rank,48) == 0:
        if rank == 0:
            file = open(pickle_file, 'wb')
            #pickle.dump((image, time_val, particle_positions, particle_masses), file)
            pickle.dump((particle_x_pos, particle_y_pos), file)
            file.close()
            print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
        del particle_x_pos
        del particle_y_pos
        gc.collect()
        
print('FINISHED GETTING PARTICLE DATA ON RANK', rank)

proj_pickles = glob.glob(save_dir + "movie_frame_*_proj.pkl")
part_pickles = glob.glob(save_dir + "movie_frame_*_part.pkl")

x = np.linspace(-1*units['length_unit']/2., units['length_unit']/2., 800)
y = np.linspace(-1*units['length_unit']/2., units['length_unit']/2., 800)
xlim = [x[0], x[-1]]
ylim = [y[0], y[-1]]
X, Y = np.meshgrid(x, y)

for pick_it in range(len(proj_pickles)):
    save_name = save_dir + "movie_frame_" + ("%06d" % pick_it)
    if os.path.exists(save_name+".jpg") == False:
        file = open(proj_pickles[pick_it], 'rb')
        image, time_val = pickle.load(file)
        file.close()
        
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_xlabel("$x$ (AU)", labelpad=-1, fontsize=args.text_font)
        ax.set_ylabel("$y$ (AU)", fontsize=args.text_font) #, labelpad=-20
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
        plt.gca().set_aspect('equal')
        cbar = plt.colorbar(plot, pad=0.0)
        
        file = open(part_pickles[pick_it], 'rb')
        particle_x_pos, particle_y_pos = pickle.load(file)
        file.close()
        
        ax.scatter(particle_x_pos, particle_y_pos, color='c', size=0.5)
        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)

        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
        
        plt.savefig(save_name + ".jpg", format='jpg', bbox_inches='tight')
        plt.savefig(save_name + ".pdf", format='pdf', bbox_inches='tight')
        print('Created frame on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
        

