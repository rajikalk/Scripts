import process_stellar as ps
import glob
import os
import astropy.io.fits as pyfits
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

S = Simbad()
S.add_votable_fields('rv_value')

templates_dir = '/Users/rajikak/tools/templates/'
files = glob.glob('/Users/rajikak/Observational_Data/RV_standards/PyWiFeS/*/*p08.fits')
image_dir = '/Users/rajikak/Observational_Data/RV_standards/Images/'
RV_standard = []
RV = []
Spec_type = []
RA = []
DEC = []

f = open('/Users/rajikak/Observational_Data/RV_standards/RV_standard.csv', 'w')
f.write('Name,RV(km/s),Spectral Type,RA,DEC\n')

for file in files:
    print "Doing file:", file
    hdu = pyfits.open(file)
    Obj_name = hdu[0].header['OBJNAME']
    if Obj_name == 'object' or Obj_name == '':
        Obj_name = hdu[0].header['OBJECT']
    if not os.path.exists(image_dir+Obj_name +'/'+file.split('/')[-1]+'.png'):
        if not os.path.exists(image_dir + Obj_name):
            os.makedirs(image_dir + Obj_name)

        flux,wave = ps.read_and_find_star_p08(file, fig_fn=image_dir+Obj_name +'/'+file.split('/')[-1]+'.png')
        
        if Obj_name not in RV_standard:
            RV_standard.append(Obj_name)
            
            custom_simbad = Simbad()
            custom_simbad.add_votable_fields("sptype")
            custom_simbad.add_votable_fields("rv_value")
            data = custom_simbad.query_object(Obj_name)
            
            try:
                result = Vizier.query_object(Obj_name)
                interesting_table = result['J/ApJS/141/503/table1']
                object_RV = float(interesting_table[0]['__RV_'])/1000.
            except:
                object_RV = float(data['RV_VALUE'])
            
            sptype = data['SP_TYPE'][-1]
            Spec_type.append(sptype)
            RV.append(object_RV)
            ra = str(data['RA'][-1])
            RA.append(ra)
            dec = str(data['DEC'][-1])
            DEC.append(dec)
            write_string = Obj_name + ',' + str(object_RV) + ',' + sptype + ',' + ra + ',' + dec + '\n'
            f.write(write_string)
        else:
            ind = RV_standard.index(Obj_name)
            object_RV = RV[ind]
            
        print "Object:", Obj_name, "has RV:", object_RV

        ps.make_wifes_p08_template(file, templates_dir, rv=object_RV)

f.close()