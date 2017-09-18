import os
import glob
from astropy.io import fits as pyfits

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sin", "--single_arc", help="Do you want to use only the first arc of the night? or two arcs on either side?", default='True')
    args = parser.parse_args()
    return args

args = parse_inputs()

files = glob.glob('T2m3w*.fits')

bias_files = []
flat_files = []
arc_files = []
object_files = []

for file in files:
    hdu = pyfits.open(file)
    imtyp = hdu[0].header['IMAGETYP']
    if imtyp == 'zero':
        bias_files.append(file.split('.fits')[0])
    elif imtyp == 'flat':
        flat_files.append(file.split('.fits')[0])
    elif imtyp == 'arc':
        arc_files.append(file.split('.fits')[0])
    elif imtyp == 'object':
        object_files.append(file.split('.fits')[0])

f = open('save_red_metadata_arcs.py', 'w')
f.write('import pickle\n')
f.write('\n')
f.write('#------------------------------------------------------\n')
f.write("bias_obs = ['" + "',\n            '".join(bias_files) + "']\n")
f.write('\n')
f.write("domeflat_obs = ['" + "',\n                '".join(flat_files) + "']\n")
f.write('\n')
f.write('twiflat_obs = []\n')
f.write('\n')
f.write('dark_obs = []\n')
f.write('\n')
if args.single_arc == 'False':
    f.write("arc_obs =  ['" + "',\n            '".join(arc_files) + "']\n")
else:
    f.write("arc_obs =  ['" + arc_files[0] + "']\n")
f.write('\n')
f.write('wire_obs = []\n')
f.write('\n')
f.write('#------------------------------------------------------\n')
f.write('sci_obs = [\n')
for arc_no in range(len(arc_files)):
    if args.single_arc == 'False':
        if arc_no == 0:
            f.write("    #"+arc_files[arc_no]+"\n")
            f.write("    {'sci'  : ['" + arc_files[arc_no] + "'],\n")
            f.write("     'arc'  : ['" + arc_files[1] + "'],\n")
            f.write("     'flat' : [],\n")
            f.write("     'sky'  : [],\n")
            f.write("     'wire' : []},\n")
        elif arc_files[arc_no] == arc_files[-1]:
            f.write("    #"+arc_files[arc_no]+"\n")
            f.write("    {'sci'  : ['" + arc_files[arc_no] + "'],\n")
            f.write("     'arc'  : ['" + arc_files[-2] + "'],\n")
            f.write("     'flat' : [],\n")
            f.write("     'sky'  : [],\n")
            f.write("     'wire' : []},\n")
        else:
            f.write("    #"+arc_files[arc_no]+"\n")
            f.write("    {'sci'  : ['" + arc_files[arc_no] + "'],\n")
            f.write("     'arc'  : ['" + arc_files[arc_no-1] + "',\n")
            f.write("               '" + arc_files[arc_no+1] + "'],\n")
            f.write("     'flat' : [],\n")
            f.write("     'sky'  : [],\n")
            f.write("     'wire' : []},\n")
    else:
        f.write("    #"+arc_files[arc_no]+"\n")
        f.write("    {'sci'  : ['" + arc_files[arc_no] + "'],\n")
        f.write("     'arc'  : ['" + arc_files[0] + "'],\n")
        f.write("     'flat' : [],\n")
        f.write("     'sky'  : [],\n")
        f.write("     'wire' : []},\n")
f.write("    ]\n")
f.write('\n')
f.write('#------------------------------------------------------\n')
f.write('std_obs = [\n')
f.write('    ]\n')
f.write('\n')
f.write('#------------------------------------------------------\n')
f.write('night_data = {\n')
f.write("    'bias' : bias_obs,\n")
f.write("    'domeflat' : domeflat_obs,\n")
f.write("    'twiflat' : twiflat_obs,\n")
f.write("    'dark' : dark_obs,\n")
f.write("    'wire' : wire_obs,\n")
f.write("    'arc'  : arc_obs,\n")
f.write("    'sci'  : sci_obs,\n")
f.write("    'std'  : std_obs}\n")
f.write('\n')
f.write("f1 = open('wifesR_re_arc_metadata.pkl', 'w')\n")
f.write("pickle.dump(night_data, f1)\n")
f.write("f1.close()\n")
f.write('\n')
f.close()