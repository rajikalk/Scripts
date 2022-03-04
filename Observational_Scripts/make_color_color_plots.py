import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", "--input_file", help="what file do you want to read in?", type=str, default='/Users/rajikak/Observational_Data/color_color_plots.csv')
    parser.add_argument("-obs", "--obs_data", help="which file has all the data with the H-alpha widths?", type=str, default='/Users/rajikak/Observational_Data/WiFeS_target_spreadsheet_H_alpha.csv')
    parser.add_argument("-save_file", "--save_file_dir", help="What do you want to save the figure as?", type=str, default='/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/color_color_plot.eps')
    args = parser.parse_args()
    return args

args = parse_inputs()

#set up lists and arrays
Object = [[],[]]
Disk_tag = [[],[]]
SpT = [[],[]]
SpT_err = [[],[]]
H_alpha = [[],[]]
H_alpha_variation = [[],[]]
J_mag = [[],[]]
J_mag_err = [[],[]]
K_mag = [[],[]]
K_mag_err = [[],[]]
W2_mag = [[],[]]
W2_mag_err = [[],[]]
W3_mag = [[],[]]
W3_mag_err = [[],[]]
W4_mag = [[],[]]
W4_mag_err = [[],[]]
K_W3 = [[],[]]
K_W3_err = [[],[]]
K_W4 = [[],[]]
K_W4_err = [[],[]]


#read in current data
print("Reading in current spreadsheet")
header = 0
with open(args.input_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if row[0] != 'UCAC4-1089654199':
                if row[1] == 'US':
                    reg_ind = 0
                else:
                    reg_ind = 1
                Object[reg_ind].append(row[0])
                Disk_tag[reg_ind].append(row[41])
                SpT[reg_ind].append(float(row[30]))
                SpT_err[reg_ind].append(float(row[31]))
                K_W3[reg_ind].append(float(row[34]))
                K_W3_err[reg_ind].append(float(row[35]))
                K_W4[reg_ind].append(float(row[38]))
                K_W4_err[reg_ind].append(float(row[39]))
                J_mag[reg_ind].append(float(row[25]))
                J_mag_err[reg_ind].append(float(row[26]))
                K_mag[reg_ind].append(float(row[28]))
                K_mag_err[reg_ind].append(float(row[29]))
                W2_mag[reg_ind].append(float(row[17]))
                W2_mag_err[reg_ind].append(float(row[18]))
                W3_mag[reg_ind].append(float(row[19]))
                W3_mag_err[reg_ind].append(float(row[20]))
                W4_mag[reg_ind].append(float(row[21]))
                W4_mag_err[reg_ind].append(float(row[22]))
        if header == 0:
            header = 1
    f.close()

H_alpha = np.array([np.zeros(len(Object[0])), np.zeros(len(Object[1]))])
H_alpha_variation = np.array([np.zeros(len(Object[0])), np.zeros(len(Object[1]))])

header = 0
with open(args.obs_data, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if row[0] != 'UCAC4-1089654199':
                if row[1] == 'US':
                    reg_ind = 0
                else:
                    reg_ind = 1
                if row[0] in Object[reg_ind]:
                    ind = Object[reg_ind].index(row[0])
                    H_alpha[reg_ind][ind] = float(row[7])
                    H_alpha_variation[reg_ind][ind] = float(row[8])
        if header == 0:
            header = 1
    f.close()

xlabels = ['G0','G2','G5','G8','K0','K2','K4','K5','K7','M0','M1','M2','M3','M4','M5','M6','M7']
J_K_color = [0.33,0.34,0.38,0.45,0.49,0.56,0.65,0.69,0.76,0.86,0.86,0.87,0.87,0.88,0.92,0.99,1.07]

J_mag = np.array([np.array(J_mag[0]),np.array(J_mag[1])])
J_mag_err = np.array([np.array(J_mag_err[0]),np.array(J_mag_err[1])])
K_mag = np.array([np.array(K_mag[0]),np.array(K_mag[1])])
K_mag_err = np.array([np.array(K_mag_err[0]),np.array(K_mag_err[1])])
W3_mag = np.array([np.array(W3_mag[0]),np.array(W3_mag[1])])
W3_mag_err = np.array([np.array(W3_mag_err[0]),np.array(W3_mag_err[1])])
W4_mag = np.array([np.array(W4_mag[0]),np.array(W4_mag[1])])
W4_mag_err = np.array([np.array(W4_mag_err[0]),np.array(W4_mag_err[1])])
J_K_min = 0.3
J_K_max = 1.3
J_K = np.arange(J_K_min, J_K_max, 0.01)
J_K_W2_1 = np.arange(J_K_min, 0.86, 0.01)
J_K_W2_2 = np.arange(0.86, J_K_max, 0.01)

plt.clf()
fig, (ax) = plt.subplots(1, 1)
fig.set_size_inches(6,5)
plt.axhline(y=0.0, color='k', ls='--')
plt.errorbar(SpT[0], H_alpha[0], xerr=SpT_err[0], yerr=H_alpha_variation[0], fmt='bo', label='Upper Scorpius')
plt.errorbar(SpT[1], H_alpha[1], xerr=SpT_err[1], yerr=H_alpha_variation[1], fmt='r^', label='Upper Centaurus-Lupus')
#plt.scatter((J_mag[0]-K_mag[0]), H_alpha[0], color='b', marker='o', label='Upper Scorpius')
#plt.scatter((J_mag[1]-K_mag[1]), H_alpha[1], color='r', marker='^', label='Upper Centaurus-Lupus')
plt.ylabel('EW(H'+r'$ \alpha $'+')')
plt.xlabel('J-K')
#plt.ylim([-45.,5.])
#ax.set_xticks(J_K_color)
#ax.set_xticklabels(xlabels)
#ax.xaxis.set_major_locator(FixedLocator(J_K_color))
plt.legend(loc='lower left')
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/H_alpha_EW_variation.eps', bbox_inches='tight', pad_inches = 0.02)

plt.clf()
fig, (ax) = plt.subplots(1, 1)
fig.set_size_inches(6,5)
plt.axhline(y=0.0, color='k', ls='--')
plt.errorbar(K_W3[0], H_alpha[0], xerr=K_W3_err[0], yerr=H_alpha_variation[0], fmt='bo', label='Upper Scorpius')
plt.errorbar(K_W3[1], H_alpha[1], xerr=K_W3_err[1], yerr=H_alpha_variation[1], fmt='r^', label='Upper Centaurus-Lupus')
#plt.scatter((J_mag[0]-K_mag[0]), H_alpha[0], color='b', marker='o', label='Upper Scorpius')
#plt.scatter((J_mag[1]-K_mag[1]), H_alpha[1], color='r', marker='^', label='Upper Centaurus-Lupus')
plt.ylabel('EW(H'+r'$ \alpha $'+')')
plt.xlabel('K-W3')
#plt.ylim([-45.,5.])
#ax.set_xticks(J_K_color)
#ax.set_xticklabels(xlabels)
#ax.xaxis.set_major_locator(FixedLocator(J_K_color))
plt.legend(loc='lower left')
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/H_alpha_EW_variation_W3.eps', bbox_inches='tight', pad_inches = 0.02)

'''
plt.clf()
fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True)
fig.set_size_inches(7,8)
ax1.plot(J_K, ((1.1*J_K)-0.3),'k--')
ax1.errorbar(SpT[0], K_W3[0], xerr=SpT_err[0], yerr=K_W3_err[0], fmt='bo', label='Upper Scorpius')
ax1.errorbar(SpT[1], K_W3[1], xerr=SpT_err[1], yerr=K_W3_err[1], fmt='r^', label='Upper Centaurus-Lupus')
#ax1.scatter((J_mag[0]-K_mag[0]), (K_mag[0]-W3_mag[0]), color='b', marker='+', label='Upper Scorpius')
#ax1.scatter((J_mag[1]-K_mag[1]), (K_mag[1]-W3_mag[1]), color='r', marker='^', label='Upper Centaurus-Lupus')
ax1.set_ylabel('K-W3')
ax1.set_ylim([-0.5,4.])
ax1.set_xlim([J_K_min,J_K_max])
ax1.legend(loc='upper left')
plt.setp([ax1.get_yticklabels()[::2]], visible=False)
ax2.plot(J_K, ((2.08*J_K)-0.48),'k--')
ax2.errorbar(SpT[0], K_W4[0], xerr=SpT_err[0], yerr=K_W4_err[0], fmt='bo', label='Upper Scorpius')
ax2.errorbar(SpT[1], K_W4[1], xerr=SpT_err[1], yerr=K_W4_err[1], fmt='r^', label='Upper Centaurus-Lupus')
#ax2.scatter((J_mag[0]-K_mag[0]), (K_mag[0]-W4_mag[0]), color='b', marker='+', label='Upper Scorpius')
#ax2.scatter((J_mag[1]-K_mag[1]), (K_mag[1]-W4_mag[1]), color='r', marker='^', label='Upper Centaurus-Lupus')
ax2.set_ylabel('K-W4')
ax2.set_xlabel('J-K')
ax2.set_xlim([J_K_min,J_K_max])
plt.setp([ax2.get_yticklabels()[-1]], visible=False)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
fig.subplots_adjust(hspace=0)
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/color_color_plot_2.eps', bbox_inches='tight', pad_inches = 0.02)
'''

plt.clf()
fig, (ax1) = plt.subplots(1, 1,sharex=True)
ax1.plot(J_K, ((2.08*J_K)-0.48),'k--')
ax1.errorbar(SpT[0], K_W4[0], xerr=SpT_err[0], yerr=K_W4_err[0], fmt='bo', label='Upper Scorpius')
ax1.errorbar(SpT[1], K_W4[1], xerr=SpT_err[1], yerr=K_W4_err[1], fmt='r^', label='Upper Centaurus-Lupus')
ax1.legend(loc='upper left')
ax1.set_ylabel('K-W4')
ax1.set_xlabel('J-K')
ax1.set_xlim([J_K_min,J_K_max])
fig.subplots_adjust(hspace=0)
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/color_color_plot_2.eps', bbox_inches='tight', pad_inches = 0.02)

plt.clf()
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True)
fig.set_size_inches(6,9)
ax1.plot(J_K_W2_1, np.ones(len(J_K_W2_1))*0.21, 'k--')
ax1.plot(J_K_W2_2, ((4.17*J_K_W2_2)-3.3762), 'k--')
ax1.scatter(SpT[0], (K_mag[0]-W2_mag[0]), color='b', marker='+', label='Upper Scorpius')
ax1.scatter(SpT[1], (K_mag[1]-W2_mag[1]), color='r', marker='^',label='Upper Centaurus-Lupus')
ax1.set_ylabel('K-W2')
ax1.set_xlim([J_K_min,J_K_max])
ax2.plot(J_K, ((1.1*J_K)-0.3),'k--')
ax2.scatter(SpT[0], (K_mag[0]-W3_mag[0]), color='b', marker='+', label='Upper Scorpius')
ax2.scatter(SpT[1], (K_mag[1]-W3_mag[1]), color='r', marker='^', label='Upper Centaurus-Lupus')
ax2.set_ylabel('K-W3')
ax2.set_xlim([J_K_min,J_K_max])
ax3.plot(J_K, ((2.08*J_K)-0.48),'k--')
ax3.scatter(SpT[0], (K_mag[0]-W4_mag[0]), color='b', marker='+', label='Upper Scorpius')
ax3.scatter(SpT[1], (K_mag[1]-W4_mag[1]), color='r', marker='^', label='Upper Centaurus-Lupus')
ax3.set_ylabel('K-w4')
ax3.set_xlabel('J-K')
ax3.set_xlim([J_K_min,J_K_max])
plt.setp([ax2.get_yticklabels()[-1]], visible=False)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
plt.setp([ax3.get_yticklabels()[-1]], visible=False)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
fig.subplots_adjust(hspace=0)
plt.savefig('/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/color_color_plot_3.eps', bbox_inches='tight', pad_inches = 0.02)
