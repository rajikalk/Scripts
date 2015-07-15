import csv
import pyfits
import commands

#create input list:
fits_list = commands.getoutput("ls T2m3wr*.fits").split("\n")

#read in objects:
objects = []
with open("log_20150331.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        index_str = ("%04d" % int(row[0]))
        objects.append([index_str, row[1]])

ob_it = 0
fit_it = 0
while fit_it < len(fits_list):
    while ob_it < len(objects):
        f = pyfits.open(fits_list[fit_it], mode='update')
        if f[0].header['OBJECT'] == '':
            if fits_list[fit_it][-9:-5] == objects[ob_it][0]:
                f[0].header['OBJECT']=objects[ob_it][1]
                f.flush()
                fit_it = fit_it + 1
                print "found fit_it =", fit_it
            else:
                ob_it = ob_it + 1
        else:
            fit_it = fit_it + 1