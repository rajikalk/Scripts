import csv
import numpy as np

#read in table
Object = []
Write_lines = []
header = 0
with open('Table_2.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            Object.append(row[0] + '&' + row[1] + '&' + row[2] + '&' + row[3] + '&' + row[4])
            counter = 5
            Write_lines.append([])
            while counter < len(row):
                if row[counter] != '':
                    line = row[counter+1] + '&' + str((np.round(float(row[counter+2])*100.))/100.) + '&' + str((np.round(float(row[counter+3])*100.))/100.) + '\\'
                    Write_lines[-1].append(line)
                counter = counter + 6
        if header == 0:
            header =1
f.close()

f = open('Table_2.tex', 'w')

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
f.write('        \caption{Full list of objects observed in this survey with each epoch of data. The spectral type listed is the spectral type of the preferred template used for each object. The error on the radial velocity is the error from the cross-correlation and the template precision (see Table. 3 in Paper) added in quadrature.}\n')
f.write('        \\')
f.write('begin{tabular}{|c|c|c|c|c|c|c|c|}\n')
f.write('            Object & Region & RA & DEC & Template SpT. & MJD & $v_r$ (km/s) & $\sigma_{v_r}$ (km/s) \\')
f.write('\\')
f.write('\n')
f.write('                \\hline\n')

for obj in range(len(Object)):
    for epoch in range(len(Write_lines[obj])):
        if epoch == 0:
            f.write(Object[obj] + '&' + Write_lines[obj][epoch])
        else:
            f.write('&&&&&' + Write_lines[obj][epoch])
        f.write('\\')
        f.write('\n')
    f.write('    \\hline\n')

f.write('        \end{tabular}\n')
f.write('        \label{table:targets}\n')
f.write('\end{table*}\n')
f.write('\n')
f.write('\end{document}\n')
f.close()
