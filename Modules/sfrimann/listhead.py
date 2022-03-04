#!/usr/bin/python

import sys
import pyfits


def main(argv):

  if len(argv) != 2:
	print "Usage : listhead <filename> <extention>"
	sys.exit (1)

  try:
    f = pyfits.open(argv[0])
    hdr = f[int(argv[1])].header
    print hdr	
    f.close()
  except Exception,e:
    print "Error : " + str(e)

if __name__ == "__main__":
	main(sys.argv[1:])
