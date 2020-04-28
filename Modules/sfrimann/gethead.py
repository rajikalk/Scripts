#!/usr/bin/python

import sys
from pyfits import getval


def main(argv):

  if len(argv) < 3:
	print "Usage : gethead <filename> <keyw1> [<keyw2 keyw3 keywn>] <#extension>"
	sys.exit (1)

  try:
    file = argv[0]
    ext  = int(argv[-1])
    for key in argv[1:-1]:
    	v = getval(file,key,ext)    	
        print v,
  except Exception,e:
    print "Error : " + str(e)	
	

if __name__ == "__main__":
	main(sys.argv[1:])
