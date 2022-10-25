#Plots MESA luminsoties calculated for all sinks int eh G100 512 run

import matplotlib.pyplot as plt
import numpy as np
import mesaPlot
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
m=mesaPlot.MESA()
p=mesaPlot.plot()
lsun = 3.828e26*1e7 # solar luminosity in erg

rank = CW.Get_rank()
size = CW.Get_size()

def plot_luminosity(sink=5, max_age=150000):
    m.log_fold='/scratch/troels/IMF_512/mesa/sink_{:04d}/LOGS'.format(sink)
    if not os.path.isdir(m.log_fold):
        return
    m.loadHistory()
    mass = m.hist.star_mass
    age = m.hist.star_age
    idx = np.where(age <= max_age)
    age = age[idx]
    lum = 10.**m.hist.log_L
    lacc = m.hist.extra_lum / lsun
    ltot = lum + lacc

    fig, ax = plt.subplots(2,figsize=(15,15))

    ax[0].plot(age,ltot[idx],label='Ltot')
    ax[0].plot(age,lum[idx],label='Lstar')
    ax[0].plot(age,lacc[idx],label='Lacc')
    ax[0].legend()
    ax[0].set_ylim(bottom=1.)
    ax[0].set_yscale('log')

    ax[1].plot(age,mass[idx],label='Mass')
    ax[1].legend()
    ax[1].set_ylim(bottom=1e-2)
    ax[1].set_yscale('log')

sink_files = sorted(glob.glob('/data/scratch/troels/IMF_512/mesa/sink_*/LOGS'))
if rank == 0:
    print('sink_files:', sink_files)
max_age=150000
rit = -1
for sink_file in sink_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        try:
            m.log_fold=sink_file
            m.loadHistory()
            mass = m.hist.star_mass
            age = m.hist.star_age
            idx = np.where(age <= max_age)
            age = age[idx]
            lum = 10.**m.hist.log_L
            lacc = m.hist.extra_lum / lsun
            ltot = lum + lacc
            
            plt.clf()
            plt.semilogy(age, lum[idx], label='L$_{star}$')
            plt.semilogy(age, lacc[idx], label='L$_{acc}$')
            plt.semilogy(age, ltot[idx], label='L$_{tot}$')
            plt.legend()
            plt.xlim([0, 150000])
            plt.ylim([np.min([np.min(lum[idx]), np.min(lacc[idx])]), np.max(ltot[idx])])
            plt.xlabel('time (yr)')
            plt.ylabel('Luminosity')
            plot_name = "luminosity_" + sink_file.split("mesa/")[-1].split("/LOGS")[0]
            plt.savefig(plot_name + ".pdf", format='pdf', bbox_inches='tight')
            print("plotted", plot_name, "on rank", rank)
        except:
            print("can't read file")

print("finished plotting on rank", rank)
