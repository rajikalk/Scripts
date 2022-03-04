#!/usr/bin/env python 

#################################
#                               #
#  ALFOSC ImStack Version 0.1   #
#  Written by: Soeren Frimann   #
#  sfrimann@gmail.com           #
#                               #
#################################

from sys import exit
import time
from subprocess import Popen
from pyfits import getval
from os import getcwd, chdir, remove, path, environ, getenv, makedirs, listdir, system, unlink, kill
from os.path import abspath, isdir, join, isfile
from glob import glob
import optparse

class FileError(Exception):
  pass

############################### Command-line parser ###############################################
parser = optparse.OptionParser(description='ALFOSC ImCombine. Software for aligning and combining CCD images', version="%prog 0.1")

parser.add_option('--logindir', action='store', dest='logindir', default=getcwd(),
                  help='Pyraf login directory')

parser.add_option('--workdir', action='store', dest='workdir', default=getcwd(),
                  help='Working directory')

parser.add_option('--align', action='store_true', dest='align', default=False,
                  help="Align Images")

parser.add_option('--combine', action='store_true', dest='combine', default=False,
                  help="Combine Images")

parser.add_option('--reference', action='store', dest='reference', default='',
                  help='Image to use as reference')

parser.add_option('--objectname', action='store', dest='objectname', default='object',
                  help='Name for the combined image')

parser.add_option('--overwrite', action='store_true', dest='overwrite', default=False,
                  help="Overwrite old files")

parser.add_option('-v','--verbose', action='store_true', dest='verbose', default=False,
                  help="Show debugging information")

parser.add_option('--check_bias', action='store_true', dest='checkbias', default=False,
                  help="Check bias level of all pictures by looking at overscan")

parser.add_option('--cleanup', action='store_true', dest='cleanup', default=False,
                  help="Cleanup everything created by ImStack from workdir")

options, remainder = parser.parse_args()

############################### Set options #######################################################
#Working directories
logindir     = abspath(options.logindir)
workdir      = abspath(options.workdir)
#Processing options
do_align     = options.align
do_combine   = options.combine
do_overwrite = options.overwrite
do_verbose   = options.verbose
do_cleanup   = options.cleanup
#Image names
reference    = options.reference
objectname   = options.objectname
#filelists
filelists    = remainder

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

def green(text):
  return teffects.green + text + teffects.reset

def bold(text):
  return teffects.bold + text + teffects.reset

def underline(text):
   return teffects.underline + text + teffects.reset

############################### DS9 ###############################################################
def detectDS9():
  from os import listdir, stat, getuid
  
  for pr in [x for x in listdir('/proc') if x.isdigit()]:
    try:
      cmdline = open("/proc/%s/cmdline" % pr).read()
      if cmdline.split('\0')[0][-3:] == "ds9":
        if stat("/proc/%s/cmdline" % pr)[4] == getuid():
          ds9_pid = int(pr)
          break
    except (OSError,IndexError):
      pass
    else:
        ds9_pid = None
  return ds9_pid

sleeptime = 3.0 # INCREASE THIS IF WORKING REMOTELY
def start_ds9():
  print 'starting ds9...'
  
  global ds9_pid, ds9_proc
  if ds9_pid == None: # start-up of first instance of ds9
    ds9_proc = Popen(["/usr/local/bin/ds9",
                      #"-zoom", "to fit"])
                      "-geometry", "620x780", 
                      "-cmap", "bb"])
    ds9_pid=ds9_proc.pid
    # give ds9 some time to start up
    print 'increase decreas Sleeptime if needed'
    time.sleep(sleeptime)
  print '...done'

############################### Function for chaning directory ####################################
def cd(target):
  chdir(target)
  iraf.cd(target)

############################### load iraf settings ################################################
def load_iraf():
  iraf.images(_doprint=0,Stdout="/dev/null")
  iraf.immatch(_doprint=0,Stdout="/dev/null")
  
  iraf.unlearn("imhead")
  iraf.unlearn("display")
  iraf.unlearn("imalign")
  iraf.unlearn("imcombine")
  
  iraf.display.fill = 'yes'
  
  iraf.images.immatch.imalign.shifts        = '@shifts.shift'
  iraf.images.immatch.imalign.boxsize       = 7
  iraf.images.immatch.imalign.bigbox        = 11
  iraf.images.immatch.imalign.niterate      = 3
  iraf.images.immatch.imalign.interp_type   = 'linear'
  iraf.images.immatch.imalign.boundary_type = 'nearest'
  iraf.images.immatch.imalign.trimimages    = 'yes'
  iraf.images.immatch.imalign.verbose       = 'yes'
  
  iraf.images.immatch.imcombine.combine  = 'average'
  iraf.images.immatch.imcombine.reject   = 'avsigclip'
  iraf.images.immatch.imcombine.scale    = 'exp'
  iraf.images.immatch.imcombine.zero     = 'none'
  iraf.images.immatch.imcombine.weight   = 'exp'
  iraf.images.immatch.imcombine.expname  = 'EXPTIME'
  iraf.images.immatch.imcombine.rdnoise  = 'RDNOISE'
  iraf.images.immatch.imcombine.gain     = 'GAIN'

############################### Start Pyraf #######################################################
currentdir = getcwd()
chdir(logindir)
from pyraf import iraf
chdir(currentdir)

############################### Set iraf settings #################################################
iraf.set(stdimage = "imt2048")
iraf.set(imtype="fits")
load_iraf()

############################### Print logo and program version ####################################
print '\n'
print header(" ALFOSC")
print header("  _____            _____                _     _             ")
print header(" |_   _|          / ____|              | |   (_)            ")
print header("   | |  _ __ ___ | |     ___  _ __ ___ | |__  _ _ __   ___  ") 
print header("   | | | '_ ` _ \| |    / _ \| '_ ` _ \| '_ \| | '_ \ / _ \ ")
print header("  _| |_| | | | | | |___| (_) | | | | | | |_) | | | | |  __/ ")
print header(" |_____|_| |_| |_|\_____\___/|_| |_| |_|_.__/|_|_| |_|\___| ")

print header("\n Version: 0.1\n")

############################### Debugging #########################################################
if do_verbose:
  print "logindir      : ", logindir
  print "workdir       : ", workdir
  print "do_align      : ", do_align
  print "do_combine    : ", do_combine
  print "do_overwrite  : ", do_overwrite
  print "do_verbose    : ", do_verbose
  print "do_cleanup    : ", do_cleanup
  print "reference     : ", reference
  print "objectname    : ", objectname
  print "filelists     : ", filelists
  print ""

############################### Cleanup ###########################################################
if do_cleanup:
  print "Cleaning ", workdir
  chdir(workdir)
  system("rm alg*.fits")
  system("rm alginput.list")
  system("rm algoutput.list")
  system("rm images.coord")
  system("rm images.shift")
  system("rm combineinput.list")
  print "Cleanup completed"
  exit(0)

############################### Sanity checks #####################################################
if len(filelists) == 0:
  raise FileError("Didn't give any lists as input")

for filelist in filelists:
  if not isfile(join(workdir,filelist)):
    raise FileError("Input file %s doesn't exist") % filelist

if reference and not isfile(reference):
  raise FileError("Error: Reference file %s does not exist" % reference)

############################### Start DS9 #########################################################

ds9_pid = None
start_ds9()

############################### Do Business #######################################################
cd(workdir)

#Read input files
files = []
for filelist in filelists:
  if filelist[-5:] == '.fits':
    files.append(filelist)
  else:
    files.extend(open(filelist).read().split())

if len(files) == 0:
  raise FileError("Error: No images in filelists. Exiting")

print files

#Display images
for i,afile in enumerate(files):
  if not isfile(afile):
    raise FileError("Error: Image %s does not exist" % afile)
  iraf.display.frame = i+1
  iraf.display(afile)
  iraf.imhead(afile)

############################### Aligning ##########################################################
if do_align and len(files) > 1:
  #Determine reference image
  if not reference:
    print 'Image number  ,  Header info'
    for k,afile in enumerate(files):
      head = iraf.imhead(afile,Stdout=1)
      print str(k)+'             ', head[0]
    print green("Choose which image to use as reference, by typing in it's number.")
    print green("Note that due to a bug in iraf.imalign you should use the smallest")
    print green("image as reference.\n")
    number = raw_input("Give the number of the image you want to use as reference ")
    reference = files[int(number)]
  
  #input/output lists
  fin  = open('alginput.list','w+')
  fout = open('algoutput.list','w+')
  for afile in files:
    print >> fin, afile
    print >> fout, 'alg' + afile
  fin.close()
  fout.close()
  
  print green("\nUse imexam to mark stars in your reference image.")
  print green("Use keystroke 'a' to measure a point source, and keystroke 'q' to quit.")
  print green("The reference image is: %s\n" % reference)
  
  #Imexam
  output = iraf.imexamine(Stdout=1)

  #make coordinate file
  f = open('images.coord','w+')
  for i in range(1,len(output)):
    splt = output[i].split()
    if len(splt) == 4:
      print >> f, splt[0] + '   ' + splt[1]
  f.close()

  print green("\nUse imexam to mark a star in all images that are to be aligned, starting with your")
  print green("reference image. Mark the same star in all images.")
  print green("Use keystroke 'a' to measure a point source, and keystroke 'q' to quit.")
  print green("The reference image is: %s\n" % reference)

  output = iraf.imexamine(Stdout=1)

  #Make shift file
  f = open('images.shift','w+')
  refx = -99.
  refy = -99.
  for i in range(1,len(output)):
    splt = output[i].split()
    if len(splt) == 4:
      if refx == -99.:
        refx = float(splt[0])
        refy = float(splt[1])
        continue
      print >> f, str(refx - float(splt[0])) + '   ' + str(refy - float(splt[1]))
  f.close()

  iraf.images.immatch.imalign(input='@alginput.list',reference=reference,coords='images.coord',output='@algoutput.list',shifts='images.shift')

  files = open('algoutput.list').read().split()
  for i,afile in enumerate(files):
    iraf.display.frame = i+1
    iraf.display(afile)

elif do_align:
  print warning("Only one image. Aligning doesn't make sense. Skipping task.")

############################### Stacking ##########################################################
if do_combine and len(files) > 1:
  if do_align:
    inputlist = 'algoutput.list'
  else:
    #input list
    fin = open('combineinput.list','w+')
    for afile in files:
      print >> fin, afile
    fin.close()
    inputlist = 'combineinput.list'
  
  iraf.images.immatch.imcombine('@'+inputlist,objectname+'.fits')
  
  #Display
  iraf.display.frame = i+2
  iraf.display(objectname)

elif do_combine:
  print warning("Only one image. Combining doesn't make sense. Skipping task.")

#kill ds9
raw_input("Press ENTER to kill ds9-display")
kill(ds9_pid,9)