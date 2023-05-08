import glob

Spins_dirs = sorted(glob.glob('/home/kuruwira/fast/Protostellar_spin/Flash_2023/*'))
Levels = ['Lref_9', 'Lref_10', 'Lref_11']
Mach_vals = ['Mach_0.1', 'Mach_0.2']
Setup = 'Binary'

header_line = 'Spin & Lref & Mach_0.1 & Mach_0.2\\'

Write_lines = []
for Lref in Levels:
    for Spin_dir in Spins_dirs:
        spin_val = Spin_dir.split('/')[-1].split('_')[-1]
        table_line = spin_val + '&' + Lref
        for Mach_val in Mach_vals:
            sink_evol_file = Spin_dir + '/' + Setup + '/' + Mach_val + '/' + Lref.split('_')[-1] + '/sinks_evol.dat'
            try:
                f1 = open(sink_evol_file, "r")
                last_lines = f1.readlines()[-3:]
                f1.close()
                last_time = last_lines[-1][12:].split('  ')[1]
                star_count = 0
                for line in last_lines[::-1]:
                    if line[12:].split('  ')[1] == last_time:
                        star_count = star_count + 1
                table_line = table_line + '&' +str(star_count)
            except:
                table_line = table_line + '&' + 'N/A'
        table_line = table_line + '\\'
        Write_lines.append(table_line)

f = open(Setup+'.tex', 'w')

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
f.write('        \caption{Stars produced for the ' + Setup+ ' setup.}\n')
f.write('        \\')
f.write('begin{tabular}{|c|c|c|c|}\n')
f.write('            Spin & Lref & Mach_0.1 & Mach_0.2\\')
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
