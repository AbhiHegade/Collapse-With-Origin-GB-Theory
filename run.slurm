#!/bin/bash -l
#SBATCH --job-name=esgb_run
#SBATCH --partition=GravityTheory
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=ah30@illinois.edu
#SBATCH --time=80:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=./out.out
#SBATCH --error=./err.err

conda activate esgb
./setup_run.py
