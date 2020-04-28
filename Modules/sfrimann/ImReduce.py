#!/usr/bin/env python 

#################################
#                               #
#  ALFOSC ImReduce Version 0.1  #
#  Written by: Soeren Frimann   #
#  sfrimann@gmail.com           #
#                               #
#################################

from sys import exit
from pyfits import getval
from os  import getcwd, chdir, remove, system, listdir
from os.path import abspath, isdir, join, isfile
from glob import glob
import optparse

class FileError(Exception):
  pass

############################### Command-line parser ###############################################
parser = optparse.OptionParser(description='ALFOSC ImReduce. Software for reducing CCD images', version="%prog 0.1")

parser.add_option('--logindir', action='store', dest='logindir', default=getcwd(),
                  help='Pyraf login directory')

parser.add_option('--workdir', action='store', dest='workdir', default=getcwd(),
                  help='Working directory')

parser.add_option('--biasdir', action='store', dest='biasdir', default="bias",
                  help='Directory containing bias images')

parser.add_option('--flatdir', action='store', dest='flatdir', default="flat",
                  help='Directory containing flat images')

parser.add_option('--objectdir', action='store', dest='objectdir', default="object",
                  help='Directory containing object images')

parser.add_option('--bias', action='store_true', dest='bias', default=False,
                  help="Do bias subtraction.")

parser.add_option('--flat', action='store_true', dest='flat', default=False,
                  help="Do flat correction.")

parser.add_option('--trim', action='store_true', dest='trim', default=False,
                  help="Trim the images.")

parser.add_option('--overscan', action='store_true', dest='overscan', default=False,
                  help="Correct overscan bias level.")

parser.add_option('--overwrite', action='store_true', dest='overwrite', default=False,
                  help="Overwrite old files")

parser.add_option('-v','--verbose', action='store_true', dest='verbose', default=False,
                  help="Show debugging information")

parser.add_option('--check_bias', action='store_true', dest='checkbias', default=False,
                  help="Check bias level of all pictures by looking at overscan")

parser.add_option('--cleanup', action='store_true', dest='cleanup', default=False,
                  help="Cleanup everything created by Imreduce from workdir")


options, remainder = parser.parse_args()

############################### Set options #######################################################
#Working directories
logindir     = abspath(options.logindir)
workdir      = abspath(options.workdir)
#Image directories
biasdir      = join(workdir,options.biasdir)
flatdir      = join(workdir,options.flatdir)
objectdir    = join(workdir,options.objectdir)
objectname   = objectdir.split('/')[-1]
#Processing options
do_bias      = options.bias
do_flat      = options.flat
do_overscan  = options.overscan
do_trim      = options.trim
do_overwrite = options.overwrite
do_verbose   = options.verbose
do_checkbias = options.checkbias
do_cleanup   = options.cleanup
#Bias and trim sections
biassection  = "[2104:2145,*]"
trimsection  = "[51:2098,2:2041]"

############################### Set text-formatting ###############################################
class teffects:
  header    = '\033[1;95m'
  blue      = '\033[94m'
  green     = '\033[92m'
  warning   = '\033[93m'
  fail      = '\033[91m'
  reset     = '\033[0;0m'
  bold      = "\033[1m"
  underline = "\033[4m"

def header(text):
   return teffects.header + text + teffects.reset

def warning(text):
  return teffects.warning + 'Warning: ' + teffects.bold + text + teffects.reset

def bold(text):
  return teffects.bold + text + teffects.reset

def underline(text):
   return teffects.underline + text + teffects.reset

############################### ALFOSC filter dictionary ##########################################
alfosc_filters = {'Open'         : 'Open',
                  'U_Bes 362_60' : 'U',
                  'B_Bes 440_100': 'B',
                  'V_Bes 530_80' : 'V',
                  'R_Bes 650_130': 'R',
                  'i_int 797_157': 'i',
                  'Own 391_7'    : 'Ostlin1'}

############################### Function for finding ALFOSC filters ###############################
def find_filter(filter_dict, filters):
  filterout = ''
  for fil in filters:
    if fil not in filter_dict:
      print warning('One or more filters are not defined in alfosc_filters dictionary')
      return 'UNDEF'
    elif (fil in filter_dict) and (filter_dict[fil] != 'Open'):
      if not filterout == '':
        print warning('More than one filter in beam')
      filterout += filter_dict[fil]
  
  if filterout == '':
    return 'Open'
  
  return filterout

############################### Function for chaning directory ####################################
def cd(target):
  chdir(target)
  iraf.cd(target)

############################### load iraf settings ################################################
def load_iraf():
  iraf.imred(_doprint=0,Stdout="/dev/null")
  iraf.ccdred(_doprint=0,Stdout="/dev/null")
  
  iraf.unlearn("imstat")
  iraf.unlearn("ccdproc")
  iraf.unlearn("flatcombine")
  iraf.unlearn("zerocombine")
  
  iraf.imstat.fields = 'image,mean,stddev'
  
  iraf.noao.imred.ccdred.ccdproc.ccdtype  = ''
  iraf.noao.imred.ccdred.ccdproc.oversca  = 'no'
  iraf.noao.imred.ccdred.ccdproc.trim     = 'no'
  iraf.noao.imred.ccdred.ccdproc.zerocor  = 'no'
  iraf.noao.imred.ccdred.ccdproc.darkcor  = 'no'
  iraf.noao.imred.ccdred.ccdproc.fixpix   = 'no'
  iraf.noao.imred.ccdred.ccdproc.flatcor  = 'no'
  iraf.noao.imred.ccdred.ccdproc.illumcor = 'no'
  iraf.noao.imred.ccdred.ccdproc.fringec  = 'no'
  iraf.noao.imred.ccdred.ccdproc.readcor  = 'no'
  iraf.noao.imred.ccdred.ccdproc.scancor  = 'no'
  iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
  iraf.noao.imred.ccdred.ccdproc.interac  = 'no'
  iraf.noao.imred.ccdred.ccdproc.biassec  = biassection
  iraf.noao.imred.ccdred.ccdproc.trimsec  = trimsection
  
  iraf.noao.imred.ccdred.zerocombine.gain    = "GAIN"
  iraf.noao.imred.ccdred.zerocombine.rdnoise = "RDNOISE"
  iraf.noao.imred.ccdred.zerocombine.ccdtype = ""
  iraf.noao.imred.ccdred.zerocombine.combine = "median"
  iraf.noao.imred.ccdred.zerocombine.reject  = "avsigclip"
  iraf.noao.imred.ccdred.zerocombine.process = "no"
  iraf.noao.imred.ccdred.zerocombine.delete  = "no"
  iraf.noao.imred.ccdred.zerocombine.scale   = "none"
  
  iraf.noao.imred.ccdred.flatcombine.combine = "median"
  iraf.noao.imred.ccdred.flatcombine.reject  = "avsigclip"
  iraf.noao.imred.ccdred.flatcombine.ccdtype = ""
  iraf.noao.imred.ccdred.flatcombine.process = "no"
  iraf.noao.imred.ccdred.flatcombine.subsets = "no"
  iraf.noao.imred.ccdred.flatcombine.scale   = "mode"
  iraf.noao.imred.ccdred.flatcombine.gain    = "GAIN"
  iraf.noao.imred.ccdred.flatcombine.rdnoise = "RDNOISE"


############################### Start Pyraf #######################################################
currentdir = getcwd()
chdir(logindir)
from pyraf import iraf
chdir(currentdir)

############################### Set iraf settings #################################################
iraf.set(imtype="fits")
load_iraf()

############################### Print logo and program version ####################################
print '\n'
print header(" ALFOSC")
print header("  _____           _____          _                  ")
print header(" |_   _|         |  __ \        | |                 ")
print header("   | |  _ __ ___ | |__) |___  __| |_   _  ___  ___  ")
print header("   | | | '_ ` _ \|  _  // _ \/ _` | | | |/ __|/ _ \ ")
print header("  _| |_| | | | | | | \ \  __/ (_| | |_| | (__|  __/ ")
print header(" |_____|_| |_| |_|_|  \_\___|\__,_|\__,_|\___|\___| ")

print header("\n Version: 0.1\n")

############################### Debugging #########################################################
if do_verbose:
  print "logindir      : ", logindir
  print "workdir       : ", workdir
  print "biasdir       : ", biasdir
  print "flatdir       : ", flatdir
  print "objectdir     : ", objectdir
  print "objectname    : ", objectname
  print "do_bias       : ", do_bias
  print "do_flat       : ", do_flat
  print "do_overscan   : ", do_overscan
  print "do_trim       : ", do_trim
  print "do_overwrite  : ", do_overwrite
  print "do_verbose    : ", do_verbose
  print "do_checkbias  : ", do_checkbias
  print "do_cleanup    : ", do_cleanup
  print "biassection   : ", biassection
  print "trimsection   : ", trimsection
  print ""

############################### Cleanup ###########################################################
if do_cleanup:
  folders = listdir(workdir)
  print "Cleaning ", workdir
  for folder in folders:
    if not isdir(folder): continue
    chdir(folder)
    system("rm logfile")
    system("rm bias.fits")
    system("rm flat_*")
    system("rm rAL*.fits")
    system("rm *.list")
    system("rm *.coord")
    system("rm *.shift")
    system("rm *bAL*.fits")
    print "Cleaned up ", folder
    chdir('..')
  print "Cleanup completed"
  exit(0)

############################### Sanity checks #####################################################
if not isdir(objectdir):
  raise FileError('Object directory does not exist')

if not isdir(biasdir) and do_bias:
  raise FileError('Bias directory does not exist')

if not isdir(flatdir) and do_flat:
  raise FileError('Flat directory does not exist')

############################### Bias reduction ####################################################

if do_bias:
  cd(biasdir)

  #Name of bias file
  bias_file = "bias.fits"
  # Check if BIAS file is present. If not, make one.
  if isfile(bias_file):
    print "Bias file found : " + bold(bias_file)
    if do_overwrite:
      print "... overwriting!"
      remove(join(biasdir,bias_file))
  
  if not isfile(bias_file):
    print "Combining bias images..."
    #Do zerocombine
    iraf.noao.imred.ccdred.zerocombine("AL*.fits[1]",
                                       output = bias_file)
    
    #imstat = iraf.imstat(bias_file, Stdout=1)
    #print imstat
    print "Bias file "+bold(bias_file)+" created"
    print "... done!\n"
    
    if do_checkbias:
      print header("Overscan level of bias images")
      iraf.imstat("AL*.fits[1]"+biassection)
      print header("Overscan level of reduced bias")
      iraf.imstat(bias_file+biassection)

############################### Flat reduction ####################################################

if do_flat:
  cd(flatdir)
  
  filters = [] #Initialize filters

  filenames = glob('AL*.fits') #Find and sort filenames
  filenames.sort()
  
  #---------------------------- pre reducing ------------------------------------------------------
  if do_bias or do_overscan:
    print "Pre reducing flat images..."
    
    prereduced = glob('./bAL*.fits')
    if len(prereduced) > 0 and do_overwrite:
      print "Pre reduced images found"
      print " ... Overwriting"
    for image in prereduced:
      remove(image)
    
    if len(prereduced) == 0 or do_overwrite:
      biasboo = 'yes' if do_bias else 'no'
      overboo = 'yes' if do_overscan else 'no'
      
      fin  = open('inbias.list','w+')
      fout = open('outbias.list','w+')
      for filename in filenames:
        print >> fin, filename + '[1]'
        print >> fout, 'b' + filename
      fin.close()
      fout.close()
      
      if not isfile(join(biasdir,bias_file)) and do_bias:
        raise FileError("Couldn't find bias file: %s" % join(biasdir,bias_file))
      
      iraf.noao.imred.ccdred.ccdproc.output   = '@' + 'outbias.list'
      iraf.noao.imred.ccdred.ccdproc.overscan = overboo
      iraf.noao.imred.ccdred.ccdproc.zerocor  = biasboo
      iraf.noao.imred.ccdred.ccdproc.zero     = join(biasdir,bias_file)
      
      iraf.noao.imred.ccdred.ccdproc('@' + 'inbias.list')
    print "Pre reduction complete"
  
  #---------------------------- Check filters -----------------------------------------------------
  for filename in filenames:
    slots = [getval(filename, "ALFLTNM"),getval(filename, "FAFLTNM"),getval(filename, "FBFLTNM")] #Three filter wheels
    filters.append(find_filter(alfosc_filters,slots)) #find filters corresponding to filename

  unique_filters = list(set(filters)) #list of unique filters
  
  #---------------------------- Flat combining ----------------------------------------------------
  print "Combining flat images..."
  #Make flat for each filter
  for unique_filter in unique_filters:
    flat_file = 'flat_' + unique_filter + '.fits'
    flat_list = 'flat_' + unique_filter + '.list'

    if isfile(flat_file):
      print "Flat file found : " + bold(flat_file)
      if do_overwrite:
        print "... overwriting!"
        remove(join(flatdir,flat_file))
        remove(join(flatdir,flat_list))
    
    if not isfile(flat_file):
      filename_list = filter(lambda x: x[0] == unique_filter, zip(filters,filenames))
      
      f = open(flat_list,'w+')
      for fn in filename_list:
        if do_bias or do_overscan:
          print >> f, 'b' + fn[1]
        else:
          print >> f, fn[1] + '[1]'
      f.close()
      
      iraf.noao.imred.ccdred.flatcombine('@'+flat_list,
                                         output=flat_file)
      
      print "Flat file "+bold(flat_file)+" created"
    
    if do_checkbias:
      print header("Overscan level of reduced flats")
      iraf.imstat(flat_file+biassection)
    
  print "... done!\n"

############################### Science reduction #################################################

cd(objectdir)

reduced = glob('./rAL*.fits')
if len(reduced) > 0:
  print "Reduced images found"
  print " ... Always Overwriting"
  for image in reduced:
    remove(image)

filters = [] #Initialize filters

filenames = glob('AL*.fits') #Find and sort filenames
filenames.sort()

if do_checkbias:
  print header("Overscan level of raw images")
  iraf.imstat("AL*.fits[1]"+biassection)

for filename in filenames:
  slots = [getval(filename, "ALFLTNM"),getval(filename, "FAFLTNM"),getval(filename, "FBFLTNM")] #Three filter wheels
  filters.append(find_filter(alfosc_filters,slots)) #find filters corresponding to filename

unique_filters = list(set(filters)) #list of unique filters

biasboo = 'yes' if do_bias else 'no'
flatboo = 'yes' if do_flat else 'no'
overboo = 'yes' if do_overscan else 'no'
trimboo = 'yes' if do_trim else 'no'

if not isfile(join(biasdir,bias_file)) and do_bias:
  raise FileError("Couldn't find bias file: %s" % join(biasdir,bias_file))

for unique_filter in unique_filters:
  flat_file       = join(flatdir,"flat_" + unique_filter + ".fits")
  object_list     =  objectname + '_' + unique_filter + '.list'
  out_object_list = 'out' + objectname + '_' + unique_filter + '.list'
  
  if not isfile(flat_file) and do_flat:
    print warning("Couldn't find flat file %s. Skipping this filter" % flat_file)
    continue #starting loop from top
  
  filename_list = filter(lambda x: x[0] == unique_filter, zip(filters,filenames))
  
  fin  = open(object_list,'w+')
  fout = open(out_object_list,'w+')
  
  for fn in filename_list:
    print >> fin, fn[1] + '[1]'
    print >> fout, 'r' + fn[1]
  fin.close()
  fout.close()
  
  if do_flat:
    #Making sure no division by zero occurs
    iraf.imreplace(flat_file,value=0.00001,upper=0,radius=0)
  
  iraf.noao.imred.ccdred.ccdproc.output = '@' + out_object_list
  iraf.noao.imred.ccdred.ccdproc.overscan = overboo
  iraf.noao.imred.ccdred.ccdproc.trim     = trimboo
  iraf.noao.imred.ccdred.ccdproc.zerocor  = biasboo
  iraf.noao.imred.ccdred.ccdproc.flatcor  = flatboo
  iraf.noao.imred.ccdred.ccdproc.zero     = join(biasdir,bias_file)
  iraf.noao.imred.ccdred.ccdproc.flat     = flat_file
  
  iraf.noao.imred.ccdred.ccdproc('@' + object_list)


print bold("\nImreduce ended successfully")

###################################################################################################