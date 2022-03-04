#!/usr/bin/python
# -*- coding: utf-8 -*-

from sys import argv, exit

if len(argv) == 4:
  time1 = argv[1]
  time2 = argv[2]
  nn    = argv[3]

elif len(argv) == 3:
  time1 = argv[1]
  time2 = argv[2]
  nn    = 1.

else:
  print " Wrong number of inputs. Exiting ... \n"
  exit(1)

h1,m1,s1 = time1.split(':')
h2,m2,s2 = time2.split(':')

print ((float(h2) - float(h1))*3600 + (float(m2) - float(m1))*60 + float(s2) - float(s1))/float(nn)