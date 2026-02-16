from mpi4py.MPI import COMM_WORLD as CW
import pickle
import yt
import glob
import sys
import os
import my_ramses_fields as myf

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--out", help="Path to the folder with simulation outputs ", type=str)
    parser.add_argument("-s", "--save","--save_dir", help="Path to the folder in which to save projections ", type=str)
    
    parser.add_argument("--pv", "--projection_vectors","--projections", help="List of projection vectors (projection plane normals). ex: [1,2,3] [3,2,1] ", nargs="+")
    parser.add_argument("-c", "--sink_tag", help="Tag of the sink on which to center the projection", type=int)
    
    parser.add_argument("-w", "--width", help="Lenght of the projection in the x and y directions of the projection plane. ex: 3.4 5.2 ", type=str )
    parser.add_argument("-d", "--depth", help="Depth of the projection in the direction of the projection vector", type=list)
    parser.add_argument("-u", "--projection_unit",  help="Unit for width and depth",type=str, default="pc")
    parser.add_argument("--nv", "--north_vector", help="Specify where do you want the north of the simulation (direction of z axis) to be projected. default=[1, 0, 0] ",type=str, default=[1, 0, 0])
    
    parser.add_argument("--first", help="Specify the output from which to start projecting. ex. 'output_00187' ", type=str)
    parser.add_argument("--last", help="Specify the final output to project. ex. 'output_00187' ", type=str)
    parser.add_argument("--step", "--frequency", help="Frequency (step) of outputs on which to project. default = 1 --> makes a projection on every output from the starting and ending output.", default=1)
    
    parser.add_argument("--pp", "--parallel_projections", help="Do you want multiple projections for the same output to be done in parallel? default = false ", default="false")
    
    args = parser.parse_args()
    return args


args = parse_inputs()

#Input and save directories:
path = args.out
save_dir = args.save
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)


proj_vectors = [eval(p) for p in args.pv]
print(proj_vectors)
sink_tag = args.sink_tag
width = ((float(args.width[0]), args.projection_unit ), (float(args.width[1]),args.projection_unit))
depth = (args.depth, args.projection_unit)
projection_in_parallel = str(args.pp)
north_vector = args.nv



rank = CW.Get_rank()
size = CW.Get_size() #for parallel projections size has to be len(files)*len(projections)=n*p ; for non parallel projections size = n

#print(rank)




#all output filles
files = sorted(glob.glob(path + '/output_*'))

# if specified from_to
if (args.first):
    from_to = [args.first, args.last]
    print('fromto:',from_to)
    
    step = int(args.step)
    #finding indexes for specified filles 
    for f in range(len(files)):
        if files[f] == path+'/'+from_to[0] :
            out_min = f
        elif files[f] == path+'/'+from_to[1]:
            out_max = f 

    files = files[out_min:out_max+1:step]
    print (files , rank)


n = len(files)
p = len(proj_vectors)

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun"), "density_unit":(2998/4**3, "Msun/pc**3")}



if (projection_in_parallel == 'false') : 
    print("1 on rank", rank)
    real_file = sorted(glob.glob(files[rank] + '/info_*'))
    ds = yt.load(real_file , units_override = units_override )
    print("2")
    for p in proj_vectors:
        #Output name "/Projection_SinkTAG_#####_[#,#,#].pkl"
        pickle_file = save_dir + "/Projection_Sink" + str(sink_tag) + "_" + (str(files[rank]).split("/")[-1]) + "_" + str(p) +".pkl"

        b = None
        saved_files = sorted(glob.glob(save_dir + '/Projection_*'))

        for f in range(len(saved_files)):
            if (saved_files[f] == pickle_file):
                print(pickle_file, ' Exists!!')
                b = 1
        print("3")            
        if (b == None):
            dd=ds.all_data()
            c = [dd['sink_particle_posx'][sink_tag], dd['sink_particle_posy'][sink_tag], dd['sink_particle_posz'][sink_tag]]
            print("4")
            projection = yt.OffAxisProjectionPlot(
                ds, p, ("gas", "density"), center=c, width=width, depth=depth, north_vector=north_vector)

            with open(pickle_file, 'wb') as f:
                pickle.dump(projection, f)
        


        
        
        
        
        

elif (projection_in_parallel == 'True'):
    #Loading the same data for each projection on different ranks
    new_files=[]
    for f in range(len(files)):
        for i in range(len(proj_vectors)):
            new_files.append(files[f])
            
    #Output name "/Projection_SinkTAG_#####_[#,#,#].pkl"
    pickle_file = save_dir + "/Projection_Sink" + str(sink_tag) + "_" + (str(new_files[rank]).split("/")[-1]) + "_" + str(proj_vectors[rank%p]) +".pkl"

    saved_files = sorted(glob.glob(save_dir + '/Projection_*'))
    print("1")

    b = None

    for f in range(len(save_files)):
        if (saved_files[f] == pickle_file):
            print(pickle_file, 'Exists!!')
            b = 1 
            
    print("2")


    if (b == None):

        #loading info.txt from the correct output folder
        real_file = sorted(glob.glob(new_files[rank] + '/info_*'))

        ds = yt.load(real_file , units_override = units_override )
        

        #projecting
        dd=ds.all_data()
        c = [dd['sink_particle_posx'][sink_tag], dd['sink_particle_posy'][sink_tag], dd['sink_particle_posz'][sink_tag]]
        L = proj_vectors[rank%p]
        projection = yt.OffAxisProjectionPlot(
            ds, L, ("gas", "density"), center=c, width=width, depth=depth, north_vector=north_vector)

        image_arr = projection.frb.data[("gas", "density")]
        

        #saving
        with open(pickle_file, 'wb') as f:
            pickle.dump(projection, f)
            
