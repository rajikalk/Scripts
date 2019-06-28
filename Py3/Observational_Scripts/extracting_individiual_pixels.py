import process_stellar as ps
import sys
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import numpy as np

file = sys.argv[1]
hdu = pyfits.open(file)
h_corr = hdu[0].header['RADVEL']
Obj_name = hdu[0].header['OBJECT']
try:
    result = Vizier.query_object(Obj_name)
    interesting_table = result['J/ApJS/141/503/table1']
    object_RV = float(interesting_table[0]['__RV_'])/1000.
except:
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields("rv_value")
    data = custom_simbad.query_object(Obj_name)
    object_RV = float(data['RV_VALUE'])
x = []
y = []
y_err = []
waves = []
specs = []

for i in range(4):
    flux,wave = ps.read_and_find_star_p08(file, manual_click=True)
    spectrum,sig = ps.weighted_extract_spectrum(flux)

    templates = glob.glob('/Users/rajikak/tools/templates/*'+file.split('/')[-1].split('p08')[0]+'*')

    rv,rv_sig, temp_used = ps.calc_rv_template(spectrum,wave,sig,templates, ([0,5400],[6870,6890]), save_figures=False, save_dir='./', heliocentric_correction=h_corr)

    x.append(i)
    y.append(rv)
    y_err.append(rv_sig)


x = np.array(x)
y = np.array(y)
y_err = np.array(y_err)


plt.clf()
plt.errorbar(x, y, yerr=y_err, fmt='o')
plt.axhline(y=object_RV)
plt.axhline(y=np.median(y), linestyle='--')
plt.show()