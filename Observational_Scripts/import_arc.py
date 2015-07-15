import csv

arc_pos = []
arcs = []
with open('log_20150331.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[1] == "Ne-Ar":
            arc_pos.append(row[0])
f.close()

fout = open('save_red_metadata.py', 'r+')
with open('save_red_metadata_no_arc.py', 'r') as f:
    reader = csv.reader(f)
    writer = csv.writer(fout)
    for row in reader:
        if len(row) > 0 and row[0][0:7] == 'arc_obs':
            writer.writerow(row)
            row = reader.next()
            while row[0][0:4] == "   '":
                arcs.append(row[0][4:-1])
                writer.writerow(row)
                row=reader.next()
        if len(row) > 0 and row[0][0:7] == 'sci_obs':
            writer.writerow(row)
            row = reader.next()
            it_pos=0
            while len(row) > 0 and row[0][0:4] == "    ":
                if row[0][6:9] == 'sci':
                    it = 0
                    pos = 0
                    if row[0][-1:] == ']':
                        pos = int(row[0][-5:-2])
                    else:
                        pos = int(row[0][-4:-1])
                    if pos < 170:
                        while int(arc_pos[it])<pos:
                            it = it + 1
                    else:
                        it = len(arc_pos)-1
                    it_pos = it
                if row[0][6:9]== 'arc':
                    arc1 = arc_pos[it_pos-1]
                    arc2 = arc_pos[it_pos]
                    for i in arcs:
                        if i[-3:] == arc1:
                            arc1 = i
                        if i[-3:] == arc2:
                            arc2 = i
                    line = ['     \'arc\'  : [\'' + arc1 + '\'', '\'' + arc2 + '\']','']
                    writer.writerow(line)
                    row = reader.next()
                writer.writerow(row)
                row = reader.next()
        writer.writerow(row)
f.close()
fout.close()
