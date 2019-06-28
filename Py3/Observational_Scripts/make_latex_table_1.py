import csv
import numpy as np

#read in table
Write_lines = []
header = 0
with open('Table_1.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if float(row[10]) > 1000:
                line = row[0] + '&' + row[1] + '&' + row[2] + '&' + row[3] + '&' + row[4] + '&' + row[6] + '&' + row[7] + '&' + row[8] + '&' + str((np.round(float(row[9])*100.))/100.) + '&' + str((np.round(float(row[10])*100.))/100.) + '&' + '$>10^3$' + '(' + row[12] + ')\\'
            else:
                line = row[0] + '&' + row[1] + '&' + row[2] + '&' + row[3] + '&' + row[4] + '&' + row[6] + '&' + row[7] + '&' + row[8] + '&' + str((np.round(float(row[9])*100.))/100.) + '&' + str((np.round(float(row[10])*100.))/100.) + '&' + str((np.round(float(row[11])*100.))/100.) + '(' + row[12] + ')\\'
            Write_lines.append(line)
        if header == 0:
            header =1
f.close()

f = open('Table_1.tex', 'w')

f.write('\documentclass[useAMS,usenatbib]{mn2e}\n')
f.write('\\usepackage{cleveref}\n')
f.write('\\usepackage{natbib}\n')
f.write('\\usepackage{graphicx}\n')
f.write('\\usepackage{amsmath}\n')
f.write('\\usepackage{epstopdf}\n')
f.write('\\usepackage{subcaption}\n')
f.write('\\usepackage{float}\n')
f.write('\\')
f.write('begin{document}\n')
f.write('\n')
f.write('\\')
f.write('begin{table*}\n')
f.write('    \centering\n')
f.write('        \caption{Full list of objects observed in this survey. P(M) indicates the probability of membership to the listed region, based on Rizzuto et al. (2015) which uses proper motion and Bayesian analysis to determine a probability. ``Disk?" indicates whether the object is found to have an excess in the W3 and/or W4 WISE bands. An object has is given the label: ``YY" if it has an excess in both bands, ``YN" if it has an excess in W3 but not W4 and ``NY" if it has an excess in W4 but not W3. $N_{Obs}$ lists the number of observations obtained over the course of this survey. Temp. SpT. lists the spectral type of preferred template used to obtain radial velocities (see Section 2.3.3). $\Delta$RV is the radial velocity variation taken to be the difference between epochs with the highest and lowest radial velocity. The Bayes\' factor calculated from the Bayesian analysis described in Section 3.3. is the ratio of likelihoods of being a binary star to being a single star.}\n')
f.write('        \\')
f.write('begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
f.write('            Object & Region & RA & DEC & P(M) & Disk? & $N_{Obs}$ & Temp. SpT. & $<v_r>$ & $\Delta v_r$ & Bayes\' Factor(binary?)\\')
f.write('\\')
f.write('\n')
f.write('                \\hline\n')


for line in Write_lines:
    f.write(line)
    f.write('\\')
    f.write('\n')

f.write('    \\hline\n')
f.write('        \end{tabular}\n')
f.write('        \label{table:targets}\n')
f.write('\end{table*}\n')
f.write('\n')
f.write('\end{document}\n')
f.close()
