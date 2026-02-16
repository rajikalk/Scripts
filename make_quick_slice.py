import yt
import sys

file = sys.argv[1]
field = sys.argv[2]

if 'chk' in file:
    ds = yt.load(file)
else:
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)

slc = yt.SlicePlot(ds, 'x', field)
slc.set_width((500,'AU'))
slc.save(field+'_'+file+'.png')