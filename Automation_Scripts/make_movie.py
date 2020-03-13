#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Christoph Federrath
import os
import sys
import time
import array
import re
import subprocess
import fnmatch
from tempfile import mkstemp
from shutil import move, copyfile
from os import remove, close, makedirs, symlink
import argparse
import subprocess

# ===== MAIN Start =====

parser = argparse.ArgumentParser(description='Make movie from a list of images.')
parser.add_argument('inputfiles', type=str, nargs='+', help='image filenames')
parser.add_argument('-o', dest='outputfile', default=None, type=str, help='movie output filename')
parser.add_argument('-s', dest='step', default=1, type=int, help='step (default=1)')
parser.add_argument('ffmpeg_args', nargs=argparse.REMAINDER)

args = parser.parse_args()

# define inputfiles and outputfile
inputfiles = args.inputfiles
if args.outputfile == None:
    cur_dir = subprocess.getoutput('pwd')+'/'
    outputfile = cur_dir.split('YT_Output')[0] + 'Videos' + cur_dir.split('YT_Output')[1].split('Movie/')[0] + cur_dir.split('YT_Output')[1].split('Movie/')[1][:-1] + '.avi'
else:
    outputfile = args.outputfile

# get files in jpg tmp file list
tmpfilelist = 'jpgtmpfilelist.txt'
fileliststr = ""
for i in range(0,len(inputfiles),args.step):
	fileliststr += inputfiles[i]+" "
shellcmd = 'ls '+fileliststr+' > '+tmpfilelist
print(shellcmd)
subprocess.call(shellcmd, shell=True)

# make movie with ffmpeg
shellcmd = 'cat $(cat '+tmpfilelist+') | ffmpeg -f image2pipe -vcodec mjpeg -i - -b:v 50000k '+' '.join(args.ffmpeg_args)+' -y '+outputfile
print(shellcmd+' -pass 1')
subprocess.call(shellcmd+' -pass 1', shell=True)
print(shellcmd+' -pass 2')
subprocess.call(shellcmd+' -pass 2', shell=True)
print(outputfile+' written.')

# remove jpg tmp file list
remove(tmpfilelist)

# ===== MAIN End =====
