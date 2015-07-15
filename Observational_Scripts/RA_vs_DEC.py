import csv
import pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import itertools
import astropy.coordinates as coord
import astropy.units as u

RA = []
DEC = []
RV = []
with open("wifes_targets.txt", "rU") as f:
    reader = csv.reader(f, delimiter=' ')
    for row in reader:
        RA_val = row[1]
        DEC_val = row[2]
        RV_val = float(row[5])
        RA.append(RA_val)
        DEC.append(DEC_val)
        RV.append(RV_val)

#Conver coords to galatic coords:
for it in range(len(RA)):
    ra_val = RA[it]
    dec_val = DEC[it]
    ra = coord.Angle(ra_val, unit = u.hour)
    dec = coord.Angle(dec_val, unit = u.degree)
    RA[it] = ra.degree
    DEC[it] = dec.degree

RA_SB = RA[0:3]
DEC_SB = DEC[0:3]
RV_SB = RV[0:3]
RA = RA[4:-1]
DEC = DEC[4:-1]
RV = RV[4:-1]

#plot figure
for record, marker in itertools.izip(range(5),itertools.cycle(mlines.Line2D.filled_markers)):
    

#plot figure
fig, ax = plt.subplots()
fig.clf()
plt.scatter(RA, DEC, s=100, c=RV, cmap=cm.jet)
plt.scatter(RA_SB, DEC_SB, s=100, c=RV_SB, cmap=cm.jet, marker="^")
plt.colorbar(s, pad=0)
plt.savefig("test.png")