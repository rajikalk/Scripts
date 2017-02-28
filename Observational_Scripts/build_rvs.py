import csv
import numpy as np
import sys

new_rvs = sys.argv[1]
cur_rvs = sys.argv[2]

#define arrays: [[Object], [[MJD],[MJD]], [[RV],[RV]], [[RV_err],[RV_err]], ...]
objects = [[],[],[],[]]

#Read current RVs
counter = -1
with open(cur_rvs, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        object = row[0]
        MJD = float(row[1])
        RV = float(row[2])
        RV_err = float(row[3])
        if object != '':
            counter = counter + 1
            objects[0].append(object)
            objects[1].append([MJD])
            objects[2].append([RV])
            objects[3].append([RV_err])
        else:
            objects[1][counter].append(MJD)
            objects[2][counter].append(RV)
            objects[3][counter].append(RV_err)
f.close()

#read new RVs
with open(new_rvs, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        object = row[0]
        MJD = float(row[4])
        RV = float(row[5])
        RV_err = float(row[6])
        if object in objects[0]:
            index = objects[0].index(object)
            objects[1][index].append(MJD)
            objects[2][index].append(RV)
            objects[3][index].append(RV_err)
        else:
            objects[0].append(object)
            objects[1].append([MJD])
            objects[2].append([RV])
            objects[3].append([RV_err])
f.close()

#write new files
f = open(cur_rvs, 'r+')
for obit in range(len(objects[0])):
    for dit in range(len(objects[1][obit])):
        if dit == 0:
            f.write(objects[0][obit] + ',' + str(objects[1][obit][dit]) + ',' + str(objects[2][obit][dit]) + ',' + str(objects[3][obit][dit]) + '\n')
        else:
            f.write('' + ',' + str(objects[1][obit][dit]) + ',' + str(objects[2][obit][dit]) + ',' + str(objects[3][obit][dit]) + '\n')
f.close()