#!/usr/bin/env python
import mpi4py
import glob
import subprocess

rank = CW.Get_rank()
size = CW.Get_size()

if rank == 0:
    chk_files = sorted(glob.glob("*chk*"))
    last_chk = chk_files[-1]

    run_line = "python ~/Scripts/Automation_Scripts/prep_restart.py " + last_chk
    subprocess.run(proj_run_line, shell=True)
    subprocess.run("sb job.sh", shell=True)
