import csv
import numpy as np

#read in table
Write_lines = []
header = 0
with open('Magnitude_table.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if row[9] != '':
                line = row[0] + '&' + row[1] + '&' + row[2] + '&' + row[3] + '&' + row[5] + '&' + row[7] + '&' + row[10] + '&' + row[12] + '&' + row[14] + '\\'
            else:
                line = row[0] + '&' + row[1] + '&' + row[2] + '&' + row[3] + '&' + row[5] + '&' + row[7] + '& & & \\'
            Write_lines.append(line)
        if header == 0:
            header =1
f.close()

f = open('Mag_table.tex', 'w')

f.write('\documentclass[useAMS,usenatbib]{mn2e}\n')
f.write('\usepackage{cleveref}\n')
f.write('\usepackage{natbib}\n')
f.write('\usepackage{graphicx}\n')
f.write('\usepackage{amsmath}\n')
f.write('\usepackage{epstopdf}\n')
f.write('\usepackage{subcaption}\n')
f.write('\usepackage{float}\n')
f.write('\\')
f.write('begin{document}\n')
f.write('\n')
f.write('\\')
f.write('begin{table*}\n')
f.write('    \centering\n')
f.write('        \caption{Magnitudes of our targets. J and K magnitudes are taken from the 2MASS survey. W4 magnitude is taken from the AllWISE Survey. The G, BP and RP magnitudes are taken from the Gaia survey.}\n')
f.write('        \\')
f.write('begin{tabular}{|c|c|c|c|c|c|c|c|c|}\n')
f.write('            Object & RA & DEC & $J$ & $K$ & $W4$ & $G$ & $BP$ & $RP$')
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
