import csv
names = []
UCAC4 = []
RIK = []
TMASS = []
EPIC = []
Other = []
dates = []
#read in targets:
with open('US_Targets.csv', 'rU') as f:
    reader = csv.reader(f)
    it = 0
    for row in reader:
        if it != 0:
            names.append([row[0], row[1], row[2], row[3], row[4]])
            if len(row[0]) != 0:
                UCAC4.append(row[0])
            if len(row[1]) != 0:
                RIK.append(row[1])
            if len(row[2]) != 0:
                TMASS.append(row[2])
            if len(row[3]) != 0:
                EPIC.append(row[3])
        else:
            it = 1

print "read through observed targets"
'''
#get other names:
with open('../WiFeS_Usco/K2Campaign2targets_Sorted.csv', 'rU') as f:
    reader = csv.reader(f)
    it = 0
    for row in reader:
        if row[1] in TMASS:
            EPIC.append(row[0])
            Other.append(row[2])
            for name in names:
                if row[1] in name:
                    name.append(row[0])
                    name.append(row[2])
        else:
            for name in names:
                if len(row[2]) != 0 and row[2] in name:
                    if len(name[2]) == 0:
                        name[2] = row[1]
                        name.append(row[0])
                    else:
                        name.append(row[0])
                        EPIC.append(row[0])


print "read through K2Campaign2Targets"
'''
#go through log
with open('2015_05_log.csv', 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if len(row) > 1:
            if "/" in row[0]:
                date = row[0]
            if row[1] != "bias" or row[1] != "Flat" or row[1] != "Arc" or row[1] != "Ne-Ar" or row[1] != "Bias":
                for name in names:
                    if len(row[1]) != 0 and row[1] in name:
                        print "Object", row[1], "observed:", name
                        dates.append([name,date])

print "Read in observed targets"

#Update targets:
fout = open('US_Targets_updated.csv', 'r+')
with open('US_Targets.csv', 'rU') as f:
    reader = csv.reader(f)
    writer = csv.writer(fout)
    date_rep = 0
    for row in reader:
        line = row
        for date in dates:
            if len(row[0]) >1 and row[0] in date[0] or len(row[1]) >1 and row[1] in date[0] or len(row[2]) >1 and row[2] in date[0] or len(row[3]) >1 and row[3] in date[0] or len(row[4]) >1 and row[4] in date[0]:
                line.append(date[1])
                date_rep = date_rep + 1
                print line[-10:]
        writer.writerow(line)
'''
    for name in names:
        if len(line[0])>1 and line[0] == name[0]:
            if len(line[1]) == 0 and len(name[1]) != 0:
                line[1] = name[1]
            if len(line[2]) == 0 and len(name[3]) != 0:
                line[2] = name[3]
            if len(line[3]) == 0 and len(name[2]) != 0:
                line[3] = name[2]
            if len(line[4]) == 0 and len(name[4]) != 0:
                line[4] = name[4]
        elif len(line[1])>1 and line[1] == name[1]:
            if len(line[2]) == 0 and len(name[3]) != 0:
                line[2] = name[3]
            if len(line[3]) == 0 and len(name[2]) != 0:
                line[3] = name[2]
            if len(line[4]) == 0 and len(name[4]) != 0:
                line[4] = name[4]
'''

print "updated target list"