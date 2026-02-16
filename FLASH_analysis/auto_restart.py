#!/usr/bin/env python
from mpi4py.MPI import COMM_WORLD as CW
import glob
import subprocess

rank = CW.Get_rank()
size = CW.Get_size()

if rank == 0:
    chk_files = sorted(glob.glob("*chk*"))
    last_chk = chk_files[-1]

    run_line = "python ~/Scripts/Automation_Scripts/prep_restart.py " + last_chk
    subprocess.run(run_line, shell=True)
    subprocess.run("sbatch job.sh", shell=True)
