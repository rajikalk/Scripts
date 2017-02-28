import csv
import numpy as np
import sys
import matplotlib.pyplot as plt

J = []
K = []
W1 = []
W2 = []
W3 = []
W4 = []
J_err = []
K_err = []
W1_err = []
W2_err = []
W3_err = []
W4_err = []
Comment_flag = []

header = True
file = sys.argv[1]
band1 = sys.argv[2]
band2 = sys.argv[3]
with open(file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header == False:
            if row[7] != '':# and row[8] == 'AAAA':
                W1_mag = row[16]
                if W1_mag == '':
                    W1_mag = 0
                else:
                    W1_mag = float(W1_mag)
                W1_sigma = row[17]
                if W1_sigma == '':
                    W1_sigma = 0
                else:
                    W1_sigma = float(W1_sigma)
                W2_mag = row[18]
                if W2_mag == '':
                    W2_mag = 0
                else:
                    W2_mag = float(W2_mag)
                W2_sigma = row[19]
                if W2_sigma == '':
                    W2_sigma = 0
                else:
                    W2_sigma = float(W2_sigma)
                W3_mag = row[20]
                if W3_mag == '':
                    W3_mag = 0
                else:
                    W3_mag = float(W3_mag)
                W3_sigma = row[21]
                if W3_sigma == '':
                    W3_sigma = 0
                else:
                    W3_sigma = float(W3_sigma)
                W4_mag = row[22]
                if W4_mag == '':
                    W4_mag = 0
                else:
                    W4_mag = float(W4_mag)
                W4_sigma = row[23]
                if W4_sigma == '':
                    W4_sigma = 0
                else:
                    W4_sigma = float(W4_sigma)
                J_mag = row[24]
                if J_mag == '':
                    J_mag = 0
                else:
                    J_mag = float(J_mag)
                J_sigma = row[25]
                if J_sigma == '':
                    J_sigma = 0
                else:
                    J_sigma = float(J_sigma)
                K_mag = row[28]
                if K_mag == '':
                    K_mag = 0
                else:
                    K_mag = float(K_mag)
                K_sigma = row[29]
                if K_sigma == '':
                    K_sigma = 0
                else:
                    K_sigma = float(K_sigma)
                if 'Y' in row[7]:
                    com = 2
                elif 'M' in row[7]:
                    com = 1
                else:
                    com = 0
                J.append(J_mag)
                K.append(K_mag)
                W1.append(W1_mag)
                W2.append(W2_mag)
                W3.append(W3_mag)
                W4.append(W4_mag)
                J_err.append(J_sigma)
                K_err.append(K_sigma)
                W1_err.append(W1_sigma)
                W2_err.append(W2_sigma)
                W3_err.append(W3_sigma)
                W4_err.append(W4_sigma)
                Comment_flag.append(com)
        else:
            header = False
'''
with open(file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header == False:
            W1_mag = row[9]
            if W1_mag == '':
                W1_mag = 0
            else:
                W1_mag = float(W1_mag)
            W1_sigma = row[10]
            if W1_sigma == '':
                W1_sigma = 0
            else:
                W1_sigma = float(W1_sigma)
            W2_mag = row[11]
            if W2_mag == '':
                W2_mag = 0
            else:
                W2_mag = float(W2_mag)
            W2_sigma = row[12]
            if W2_sigma == '':
                W2_sigma = 0
            else:
                W2_sigma = float(W2_sigma)
            W3_mag = row[14]
            if W3_mag == '':
                W3_mag = 0
            else:
                W3_mag = float(W3_mag)
            W3_sigma = row[13]
            if W3_sigma == '':
                W3_sigma = 0
            else:
                W3_sigma = float(W3_sigma)
            W4_mag = row[15]
            if W4_mag == '':
                W4_mag = 0
            else:
                W4_mag = float(W4_mag)
            W4_sigma = row[16]
            if W4_sigma == '':
                W4_sigma = 0
            else:
                W4_sigma = float(W4_sigma)
            if 'NN' in row[6]:
                com = False
            else:
                com = True
            W1.append(W1_mag)
            W2.append(W2_mag)
            W3.append(W3_mag)
            W4.append(W4_mag)
            W1_err.append(W1_sigma)
            W2_err.append(W2_sigma)
            W3_err.append(W3_sigma)
            W4_err.append(W4_sigma)
            Comment_flag.append(com)
        else:
            header = False
'''
J = np.array(J)
W1 = np.array(W1)
W2 = np.array(W2)
W3 = np.array(W3)
W4 = np.array(W4)
J_err = np.array(J_err)
W1_err = np.array(W1_err)
W2_err = np.array(W2_err)
W3_err = np.array(W3_err)
W4_err = np.array(W4_err)

if band1 == 'W1':
    if band2 == 'W2':
        diff = W1 - W2
        diff_err = W1_err + W2_err
        band_arr = W2
        band_arr_err = W2_err
    if band2 == 'W3':
        diff = W1 - W3
        diff_err = W1_err + W3_err
        band_arr = W3
        band_arr_err = W3_err
    if band2 == 'W4':
        diff = W1 - W4
        diff_err = W1_err + W4_err
        band_arr = W4
        band_arr_err = W4_err
elif band1 == 'W2':
    if band2 == 'W3':
        diff = W2 - W3
        diff_err = W2_err + W3_err
        band_arr = W3
        band_arr_err = W3_err
    if band2 == 'W4':
        diff = W2 - W4
        diff_err = W2_err + W4_err
        band_arr = W4
        band_arr_err = W4_err
if band1 == 'J':
    if band2 == 'W1':
        diff = J - W1
        diff_err = J_err + W1_err
        band_arr = W1
        band_arr_err = W1_err
    if band2 == 'W2':
        diff = J - W2
        diff_err = J_err + W2_err
        band_arr = W2
        band_arr_err = W2_err
    if band2 == 'W3':
        diff = J - W3
        diff_err = J_err + W3_err
        band_arr = W3
        band_arr_err = W3_err
    if band2 == 'W4':
        diff = J - W4
        diff_err = J_err + W4_err
        band_arr = W4
        band_arr_err = W4_err
if band1 == 'K':
    if band2 == 'W1':
        diff = K - W1
        diff_err = K_err + W1_err
        band_arr = W1
        band_arr_err = W1_err
    if band2 == 'W2':
        diff = K - W2
        diff_err = K_err + W2_err
        band_arr = W2
        band_arr_err = W2_err
    if band2 == 'W3':
        diff = K - W3
        diff_err = K_err + W3_err
        band_arr = W3
        band_arr_err = W3_err
    if band2 == 'W4':
        diff = K - W4
        diff_err = K_err + W4_err
        band_arr = W4
        band_arr_err = W4_err
else:
    diff = W3 - W4
    diff_err = W3_err + W4_err
    band_arr= W4
    band_arr_err = W4_err


plotted = 0
while plotted < len(band_arr):
    for flag in range(3):
        for obj in range(len(W1)):
            if Comment_flag[obj] == flag:
                if Comment_flag[obj] == 0:
                    color = 'b'
                elif Comment_flag[obj] == 1:
                    color = 'g'
                else:
                    color = 'r'
                plt.plot((diff[obj] - diff_err[obj], diff[obj] + diff_err[obj]),(band_arr[obj], band_arr[obj]), c=color)
                plt.plot((diff[obj], diff[obj]),(band_arr[obj] - band_arr_err[obj], band_arr[obj] + band_arr_err[obj]), c=color)
                plt.plot(diff[obj], band_arr[obj], 'o', c=color)
                plotted = plotted + 1
plt.xlabel(band1 + '-' + band2 + ' (mag)')
plt.ylabel(band2 + ' (mag)')
#plt.xlim([-7,1])
#plt.ylim([-1,12])
plt.savefig(band1+'_'+band2+'.eps')
plt.clf()