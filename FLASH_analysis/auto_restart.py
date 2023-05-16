#!/usr/bin/env python
import glob
import subprocess

chk_files = sorted(glob.glob("*chk*"))
last_chk = chk_files[-1]

run_line = "python ~/Scripts/Automation_Scripts/prep_restart.py " + last_chk
subprocess.run(proj_run_line, shell=True)
subprocess.run("sb job.sh", shell=True)
