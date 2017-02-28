import process_stellar as ps
import glob
import pyfits
import sys
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

S = Simbad()
S.add_votable_fields('rv_value')

dir = sys.argv[1]
files = glob.glob(dir + '/*p08.fits')
templates = glob.glob('/home/rajikak/tools/RV_standards/*')

print "*************"

for file in files:
    a = pyfits.open(file)
    object_name = a[1].header['OBJNAME']

    #get RV data on the standard:
    #From Vizier table:
    result = Vizier.query_object(object_name)
    interesting_table = result['J/ApJS/141/503/table1']
    object_RV = float(interesting_table[0]['__RV_'])/1000.

    flux_stamp, wave = ps.read_and_find_star_p08(file)
    spectrum,sig = ps.weighted_extract_spectrum(flux_stamp)
    rv,rv_sig = ps.calc_rv_template(spectrum,wave,sig,templates, ([0,5400],[6870,6890]))
    rv += a[1].header['RADVEL']

    print "object_name: " + object_name
    print "object RV: " + str(object_RV) +"km/s"
    print "Retrieved RV from standards: "+ str(rv) + "km/s"
    print "*************"