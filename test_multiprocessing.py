import multiprocessing as mp
from multiprocessing import Process, Value, Array

def f(a, i):
    a[i] = i*i

arr = Array('i', range(10))
print arr[:]

processes = []
for i_val in range(len(arr)):
    processes.append(Process(target=f, args=(arr, i_val)))
while len(processes) > 0:
    running_processes = []
    for i in xrange(0, mp.cpu_count()):
        if len(processes) > 0:
            running_processes.append(processes.pop())
    for p in running_processes:
        p.start()
        print "starting process", processes
    for p in running_processes:
        p.join()


print arr[:]