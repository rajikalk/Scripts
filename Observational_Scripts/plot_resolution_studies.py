import csv
import matplotlib.pyplot as plt
import glob
import numpy as np

def __unicode__(self):
    return str(self.some_field) or ''

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-plot_ax", "--plot_axis", help="Do you want to plot resolution vs the no. orbs or no. systems? default='orbs'", type=str, default='orbs')
    parser.add_argument("-orbs", "--no_orbs", help="if you want to plot over the no. systems, what no. orbs do you want to use?", type=str, default='1e5')
    parser.add_argument("-sys", "--no_sys", help="if you want to plot over the no. orbs, what no. systems do you want to use?", type=str, default='1e0')
    parser.add_argument("-fig", "--fig_name", help="What do you want to save your figure as?", type=str, default='test.png')
    args = parser.parse_args()
    return args

args = parse_inputs()

if args.plot_axis == 'orbs':
    files = sorted(glob.glob(args.no_sys+'_systems/1e*/*.csv'))
else:
    files = sorted(glob.glob('1e*_systems/'+args.no_orbs+'*/*.csv'))

Objects = []
x = []
bayes_factors = []

for file in files:
    line_count = 0
    with open(file, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] not in Objects:
                Objects.append(row[0])
                bayes_factors.append([])
            obj_ind = Objects.index(row[0])
            bayes_factors[obj_ind].append(float(row[-1]))
            line_count = line_count + 1
    if line_count > 0:
        if args.plot_axis == 'orbs':
            x.append(int(float(file.split('/')[1].split('_')[0])))
        else:
            x.append(int(float(file.split('/')[0].split('_')[0])))
    print("for file:", file, "line_count =", line_count)

x = np.array(x)

bayes_factors = np.array(bayes_factors)
maxes = []
for repeat in range(len(bayes_factors[0])-1):
    maxes.append(np.max(abs(bayes_factors[:,1:] - bayes_factors[:,:-1]),axis=1))
maxes = np.array(maxes).T
norm_diff = abs(bayes_factors[:,1:] - bayes_factors[:,:-1])/maxes
means = np.mean(norm_diff,axis=0)
err = np.std(norm_diff,axis=0)

plt.clf()
for obj in range(len(Objects)):
    plt.scatter((x[1:] - x[:-1])/2., abs(bayes_factors[obj][1:] - bayes_factors[obj][:-1])/np.max(abs(bayes_factors[obj][1:] - bayes_factors[obj][:-1])))
    plt.errorbar((x[1:] - x[:-1])/2., means, yerr=err, fmt='o')
    plt.xscale('log')
plt.ylabel('Difference between iterations')
if args.plot_axis == 'orbs':
    plt.xlabel('no. of orbs')
else:
    plt.xlabel('no. of systems')
plt.savefig(args.fig_name)

