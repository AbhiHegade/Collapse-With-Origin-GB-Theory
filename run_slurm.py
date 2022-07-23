#!/usr/bin/env python

import time
from datetime import datetime
import os
import subprocess
current_time = datetime.now()
home_dir = "/home/ah30/scratch/code-f-phi"
dir_name = home_dir + "/slurm_py/" + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
with open(dir_name + "/run_py.slurm", 'w') as f:
    f.write("#!/bin/bash -l \n")
    f.write("#SBATCH --job-name=search \n")
    f.write("#SBATCH --partition=GravityTheory\n")
    f.write("#SBATCH --mail-type=FAIL\n")
    f.write("#SBATCH --mail-user=ah30@illinois.edu\n")
    f.write("#SBATCH --time=80:00:00\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks-per-node=10\n")
    f.write("#SBATCH --mem=1200MB		# memory in MB\n")
    f.write("#SBATCH --output={}/output.out\n".format(dir_name))
    f.write("#SBATCH --error={}/err.err\n".format(dir_name))
    f.write("conda activate esgb\n")
    f.write("{}/setup_run.py".format(home_dir))

subprocess.call("sbatch {}".format(dir_name + "/run_py.slurm"),shell = True)
