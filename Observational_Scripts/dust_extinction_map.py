#import numpy as np
#import h5py
import healpy as hp
import matplotlib.pyplot as plt
import pickle
import csv
from astropy import units as u
from astropy.coordinates import SkyCoord

Objects = []
Bayes_factor = []

header = 0
with open('/short/ek9/rlk100/Bayes_analysis/Mean_velocity/3_sig/bayes_factors_reviewer_response.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0] != 'UCAC4-1253626396':
            Objects.append(row[0])
            Bayes_factor.append(float(row[2]))

#set up lists and arrays
RA_US = [[],[]]
DEC_US = [[],[]]
RA_UCL = [[],[]]
DEC_UCL = [[],[]]

header = 0
with open('/home/100/rlk100/WiFeS_target_spreadsheet_reviewer_response.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            c = SkyCoord(row[2]+' '+ row[3], unit=(u.hourangle, u.deg))
            if row[0] in Objects:
                ind = Objects.index(row[0])
                if row[1] == 'US':
                    if Bayes_factor[ind] > 300.0:
                        RA_US[0].append(c.ra.value - 180)
                        DEC_US[0].append(c.dec.value)
                    else:
                        RA_US[1].append(c.ra.value - 180)
                        DEC_US[1].append(c.dec.value)
                else:
                    if Bayes_factor[ind] > 300.0:
                        RA_UCL[0].append(c.ra.value - 180)
                        DEC_UCL[0].append(c.dec.value)
                    else:
                        RA_UCL[1].append(c.ra.value - 180)
                        DEC_UCL[1].append(c.dec.value)
        if header == 0:
            header = 1

'''
def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
'''
#Lupus_x_bounds = [5, 25]
#Lupus_y_bounds = [334-180, 352-180]

#Sco_cen_y_bounds = [-58.6978, -15.6747]
#Sco_cen_x_bounds = [24.0971, 74.0141] #[204.0971-180., 254.0141-180.]


# Open the file and extract pixel information and median reddening in the far limit
'''
print "Reading in pix_info"
temp_file = open('pix_info.pkl', 'rb')
pix_info = pickle.load(temp_file)
temp_file.close()
f = h5py.File('dust-map-3d.h5', 'r')
pix_info = f['/pixel_info'][:]
temp_file = open('pix_info.pkl', 'w')
pickle.dump((pix_info), temp_file)
temp_file.close()
f.close()
print "Reading in EBV_far_median"
temp_file = open('EBV_far_median.pkl', 'rb')
EBV_far_median = pickle.load(temp_file)
temp_file.close()
print "Done reading in EBV_far_median"
print "Calculating EBV_far_median"
EBV_far_median = np.median(f['/samples'][:,:,-1], axis=1)
f.close()
temp_file = open('EBV_far_median.pkl', 'w')
pickle.dump((EBV_far_median), temp_file)
print "SAVED EBV_far_median.pkl"
temp_file.close()

# Construct an empty map at the highest HEALPix resolution present in the map
nside_max = np.max(pix_info['nside'])
n_pix = hp.pixelfunc.nside2npix(nside_max)
pix_val = np.empty(n_pix, dtype='f8')
pix_val[:] = np.nan

for nside in np.unique(pix_info['nside']):
    # Get indices of all pixels at current nside level
    idx = pix_info['nside'] == nside

    # Extract E(B-V) of each selected pixel
    pix_val_n = EBV_far_median[idx]

    # Determine nested index of each selected pixel in upsampled map
    mult_factor = (nside_max/nside)**2
    pix_idx_n = pix_info['healpix_index'][idx] * mult_factor

    # Write the selected pixels into the upsampled map
    for offset in range(mult_factor):
        pix_val[pix_idx_n+offset] = pix_val_n[:]

print "writing out pixel values"
temp_file = open('pix_vals.pkl', 'w')
pickle.dump((pix_val), temp_file)
temp_file.close()
print "wrote out pixel values"
'''
print "reading in pixel values"
temp_file = open('pix_vals.pkl', 'rb')
pix_val = pickle.load(temp_file)
temp_file.close()

lonra=[23, 75]
latra=[-61, -15]
hp.visufunc.cartview(pix_val, nest=True, xsize=4000, format=r'$%g$', title=r'$\mathrm{E} ( B-V \, )$', unit='$\mathrm{mags}$',lonra=lonra, latra=latra)
plt.savefig('test_healpy.png')
test=hp.cartview(pix_val, return_projected_map=True, lonra=lonra, latra=latra, nest=True)
plt.clf()
plt.imshow(test, origin='lower',extent=(lonra[1],lonra[0],latra[0],latra[1]), interpolation = 'none', cmap='Greys', clim=(0.0, 0.50))
plt.colorbar(pad=0.0)
plt.scatter(RA_US[0], DEC_US[0], color='b', marker='v')
plt.scatter(RA_US[1], DEC_US[1], color='b', marker='+', label='Upper Scorpius')
plt.scatter(RA_UCL[0], DEC_UCL[0], color='r', marker='v')
plt.scatter(RA_UCL[1], DEC_UCL[1], color='r', marker='+', label='Upper Centaurus-Lupus')
plt.legend(loc='lower right')
plt.xlim(lonra)
plt.ylim(latra)
plt.xlabel('RA')
plt.ylabel('DEC')
plt.savefig('test_matplotlib.eps', format='eps', bbox_inches='tight')

'''
x_pos = []
y_pos = []

for pix in pix_info:
    x_pos.append(IndexToDeclRa(int(pix['nside']),int(pix['healpix_index']))[0])
    y_pos.append(IndexToDeclRa(int(pix['nside']),int(pix['healpix_index']))[1])

x_pos = np.array(x_pos)
y_pos = np.array(y_pos)

lupus_inds = np.where((x_pos>5.)&(x_pos<25.)&(y_pos>334.)&(y_pos<352.))[0]
'''

