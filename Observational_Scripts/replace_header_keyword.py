#!/usr/bin/env python

#python script to replace a header keyword with something else.

from astropy.io import fits as pyfits

filename=raw_input('Which file needs a header replacement? :')
keyword=str.upper(raw_input('Which keyword (case insensitive)? :'))

try:
    f=pyfits.open(filename,mode='update',ignore_missing_end=True)
except Exception:
    print 'That file does not seem to exist'
    exit()

new_key=raw_input('current header keyword for '+keyword+' is set to: '+f[0].header[keyword]+'. What do you want to change it to? :')
f[0].header[keyword]=new_key

print 'Done, thanks'
f.flush()
