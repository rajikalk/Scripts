clfind2d
========

This utility provides a Python port of the clump finding algorithm developed by Jonathan Williams (see http://www.ifa.hawaii.edu/users/jpw/clumpfind.shtml and Willieam et al. 1994 for detailed documentation).


Dependences
-----------

clfind2d.py requires Scipy and Numpy to work. These modules are present by default in most Python installations. The code makes heavy use of scipy.ndimage.label function to substitute the search2d.pro (http://idlastro.gsfc.nasa.gov/) routine.


Usage
-----

To run the code on a FITS image, use

python clfind2d.py image.fits levels

where levels are the contour levels to be used (comma-separated). For example

python clfind2d.py NGC0925_FUV_NGS.fits 0.0044,0.0067,0.0089,0.0111,0.0133,
0.0156,0.0178,0.02,0.0222,0.0245,0.0267,0.0289,0.0311,0.0334,0.0356,0.0378,
0.04,0.0423,0.0445,0.0467,0.0489,0.0512,0.0534,0.0556,0.0578,0.0601,0.0623,
0.0645,0.0667,0.0689,0.0712,0.0734,0.0756,0.0778,0.0801,0.0823,0.0845,0.0867,
0.089,0.0912,0.0934,0.0956,0.0979,0.1001,0.1023,0.1045,0.1068,0.109,0.1112,
0.1134,0.1157,0.1179,0.1201,0.1223,0.1246,0.1268,0.129,0.1312,0.1335,0.1357,
0.1379,0.1401,0.1423,0.1446,0.1468,0.149,0.1512,0.1535,0.1557,0.1579,0.1601,
0.1624,0.1646,0.1668,0.169,0.1713,0.1735,0.1757,0.1779,0.1802,0.1824,0.1846,
0.1868,0.1891,0.1913,0.1935,0.1957

Another (usually more convenient) possibility is to use the the -i option to provide an interval in the shape [start,end,step]. The following command has the same effect as the previous one

python clfind2d.py -i 0.0044,0.1957,0.0023 NGC0925_FUV_NGS.fits

Use clfind2d.py -h for a complete list of all the commands available.


Performance
-----------

Although no extensive benchmark testing have been performed, the code seems to be comparable in speed and requirements to the IDL code. For levels where a high number of regions are produced, the IDL code seems to be slightly more efficient. On the contrary, runs that require the definition of a small/medium number of regions are more efficient using the Python code.

