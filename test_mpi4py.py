#!/usr/bin/env python
import mpi4py
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import glob
import pickle
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

rank = CW.Get_rank()
size = CW.Get_size()

print("Number of processes =", size)
print("Rank =", rank)

CW.Barrier()
print("-------------------------------------")
CW.Barrier()

if rank == 0:
    print("Number of processes =", size)
print("Rank =", rank)

#The confusing thing when writing scripts for parallelisation is that you need to be wear of what each individual rank sees. This can be confusing when playing with the same chuck of data.

test_array = np.zeros([4,5])
print("Rank", rank, ":")
print(test_array)

CW.Barrier()

for rank_id in range(size):
    if rank == rank_id:
        test_array[rank_id] = test_array[rank_id]+(rank_id+1)
    
CW.Barrier()

print("Rank", rank, ":")
print(test_array)

#If tasks are independant then it is easy to parallelise. You can send one independent task to each process, eg, calculate a for individual data sets?

pickle_files = glob.glob("test_data_*.pkl")

#read data:
file = open(pickle_files[rank], 'rb')
x,y = pickle.load(file)
file.close()

guess_mean = np.mean(y)
guess_std = 3*np.std(y)/(2**0.5)/(2**0.5)
guess_phase = 0
guess_freq = abs(x[np.argsort(y)[10]] - x[np.argsort(y)[0]])
guess_amp = 1
data_first_guess = guess_std*np.sin(x+guess_phase) + guess_mean
optimize_func = lambda z: z[0]*np.sin(z[1]*x+z[2]) + z[3] - y
est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
data_fit = est_amp*np.sin(est_freq*x+est_phase) + est_mean
fine_t = np.arange(0,max(x),0.1)
data_fit=est_amp*np.sin(est_freq*fine_t+est_phase)+est_mean

plt.clf()
plt.plot(x, y, '.')
plt.plot(x, data_first_guess, label='first guess')
plt.plot(fine_t, data_fit, label='after fitting')
plt.legend()
plt.savefig("fit_for_rank_"+str(rank)+".png")

#But if you can also communicate between processes. The simplest methods is to send and receive between specific ranks

data = {}

if rank == 0:
    data = {'a': 7, 'b': 3.14}
    CW.send(data, dest=1, tag=11)
elif rank == 1:
    data = CW.recv(source=0, tag=11)

print("On rank", rank, "data =", data)

CW.Barrier()

#But if your sending Numpy arrays, this method is much faster!
data = np.array([])
if rank == 0:
    data = np.arange(100, dtype=np.float64)
    CW.Send(data, dest=1, tag=13)
elif rank == 1:
    data = np.empty(100, dtype=np.float64)
    CW.Recv(data, source=0, tag=13)
    
print("On rank", rank, "data =", data)

CW.Barrier()

#But what if you want to send someone from one rank to everyone else? Well you can braodcast:

if rank == 0:
    data = {'key1' : [7, 2.72, 2+3j],
            'key2' : ( 'abc', 'xyz')}
else:
    data = None
data = CW.bcast(data, root=0)
# or CW.Bcast(data, root=0) for numpy arrays

print("On rank", rank, "data =", data)

CW.Barrier()

#You might have to do some testing of doing a task on one process and sending that to everyone else though. If you are sending big stuff it might take a while.
#Something that is a litle bit different to broadcasting but still shares information on one node with others is scatter, which can scatter elements of a list, or array to all the process. This can be good with those independant tasks. So lets try scattering the pickle files:

if rank == 0:
    pickle_file = glob.glob("test_data_*.pkl")
else:
    pickle_file = None
pickle_file = CW.scatter(pickle_file, root=0)

print("On rank", rank, "pickle_file =", pickle_file)

CW.Barrier()

#Now that each rank has a file, it can do it's own thing.

#If individual processes were doing their own tasks, but you need to collate everything, you can do that using gather, which puts everything from the different ranks in a list:

data = (rank+1)**2
print("On rank", rank, "data =", data)
CW.Barrier()

data = CW.gather(data, root=0)
print("On rank", rank, "data =", data)

CW.Barrier()

#When dealing with reading and writing files, it is probbaly best to only let one rank do this because you may run into I/O issues.
