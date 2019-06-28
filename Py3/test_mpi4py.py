

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import yt
yt.enable_parallelism()
import my_fields as myf
import sys
import time

file = sys.argv[1]
part_file = part_file = file[:-12] + 'part' + file[-5:]
ds = yt.load(file, particle_filename=part_file)

dd = ds.all_data()
start = time.time()
ds.derived_field_list
dd['Gravitational_Force_on_particles_Rad']
print('Finished deriving field list after', time.time()-start)

