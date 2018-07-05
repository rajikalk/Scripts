import numpy as np
import csv
import matplotlib.pyplot as plt
'''
import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1

mpl.rcParams['lines.linewidth']   =3
mpl.rcParams['axes.linewidth']    = 2
mpl.rcParams['xtick.major.width'] =2
mpl.rcParams['ytick.major.width'] =2
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['axes.labelsize'] = 15
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['axes.labelweight']='semibold'
mpl.rcParams['mathtext.fontset']='stix'
mpl.rcParams['font.weight'] = 'semibold'
'''
def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-mf", "--member_file", help="what file with the members do you want to read in?", type=str, default='/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/Gaia_known_members.csv')
    parser.add_argument("-file", "--input_file", help="which file has all the data from our survey?", type=str, default='/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/Gaia_our_members.csv')
    parser.add_argument("-save_file", "--save_file_dir", help="What do you want to save the figure as?", type=str, default='/Users/rajikak/Dropbox/Reggie_PhD/Papers/Multiplicity_US_UCL_2018/Resubmission_1/CMD.eps')
    args = parser.parse_args()
    return args

args = parse_inputs()

RA = [[],[]]
DEC = [[],[]]
Par = [[],[]]
Par_err = [[],[]]
G_mag = [[],[]]
G_mag_err = [[],[]]
BP_mag = [[],[]]
BP_mag_err = [[],[]]
RP_mag = [[],[]]
RP_mag_err = [[],[]]
BP_RP = [[],[]]
BP_RP_err = [[],[]]
G_mag_abs = [[],[]]
G_mag_err_abs = [[],[]]
BP_mag_abs = [[],[]]
BP_mag_err_abs = [[],[]]
RP_mag_abs = [[],[]]
RP_mag_err_abs = [[],[]]
BP_RP = [[],[]]
BP_RP_err = [[],[]]


print "reading in Known ScoCen members"
header = 0
with open(args.member_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if row[5]!= '' and row[6]!= '' and row[7] != '' and row[9] != '' and row[11] != '':
                RA[0].append(float(row[1]))
                DEC[0].append(float(row[3]))
                Par[0].append(float(row[5]))
                Par_err[0].append(float(row[6]))
                G_mag_err[0].append(1./float(row[7]))
                G_mag[0].append(float(row[8]))
                BP_mag_err[0].append(1./float(row[9]))
                BP_mag[0].append(float(row[10]))
                RP_mag_err[0].append(1./float(row[11]))
                RP_mag[0].append(float(row[12]))
        else:
            header = 1
f.close()

print "reading in our targets"
header = 0
with open(args.input_file, 'rU') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            if row[7] != '' and row[9] != '' and row[11] != '':
                RA[1].append(float(row[1]))
                DEC[1].append(float(row[3]))
                Par[1].append(float(row[5]))
                Par_err[1].append(float(row[6]))
                G_mag_err[1].append(1./float(row[7]))
                G_mag[1].append(float(row[8]))
                BP_mag_err[1].append(1./float(row[9]))
                BP_mag[1].append(float(row[10]))
                RP_mag_err[1].append(1./float(row[11]))
                RP_mag[1].append(float(row[12]))
        else:
            header = 1
f.close()

G_mag_abs[0] = np.array(G_mag[0]) + 5*(np.log10(np.array(Par[0])/1000.)+1)
G_mag_abs[1] = np.array(G_mag[1]) + 5*(np.log10(np.array(Par[1])/1000.)+1)
BP_mag_abs[0] = np.array(BP_mag[0]) + 5*(np.log10(np.array(Par[0])/1000.)+1)
BP_mag_abs[1] = np.array(BP_mag[1]) + 5*(np.log10(np.array(Par[1])/1000.)+1)
RP_mag_abs[0] = np.array(RP_mag[0]) + 5*(np.log10(np.array(Par[0])/1000.)+1)
RP_mag_abs[1] = np.array(RP_mag[1]) + 5*(np.log10(np.array(Par[1])/1000.)+1)


BP_RP[0] = np.array(BP_mag[0]) - np.array(RP_mag[0])
BP_RP_err[0] = np.sqrt(np.array(BP_mag_err[0])**2. + np.array(RP_mag_err[0])**2.)
BP_RP[1] = np.array(BP_mag[1]) - np.array(RP_mag[1])
BP_RP_err[1] = np.sqrt(np.array(BP_mag_err[1])**2. + np.array(RP_mag_err[1])**2.)

'''
BP_RP[0] = np.array(BP_mag_abs[0]) - np.array(RP_mag_abs[0])
BP_RP_err[0] = np.sqrt(np.array(BP_mag_err[0])**2. + np.array(RP_mag_err[0])**2.)
BP_RP[1] = np.array(BP_mag_abs[1]) - np.array(RP_mag_abs[1])
BP_RP_err[1] = np.sqrt(np.array(BP_mag_err[1])**2. + np.array(RP_mag_err[1])**2.)
'''
plt.clf()
#plt.scatter(BP_mag[0], BP_RP_members, color='b')
#plt.scatter(BP_mag[1], BP_RP_this_work, color='r')
plt.errorbar(BP_RP[0], G_mag_abs[0], xerr=BP_RP_err[0], yerr=G_mag_err[0], fmt='bo', label='Known Members')
plt.errorbar(BP_RP[1], G_mag_abs[1], xerr=BP_RP_err[1], yerr=G_mag_err[1], fmt='r^', label='This Work')
plt.legend(loc='best')
plt.ylabel('G$_{abs}$ (mag)')
plt.xlabel('BP - RP (mag)')
plt.gca().invert_yaxis()
plt.savefig('BP_RP_CMD.eps', bbox_inches='tight', pad_inches = 0.02)

