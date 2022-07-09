#!/usr/bin/env python
# coding: utf-8
from datetime import datetime
import os, sys, time, shutil, subprocess

class Sim:
    def __init__(self):
        self.home_dir = str(os.getcwd())
#===============================================================================
    def make_output_dir(self):
        current_time = datetime.now()
        self.output_dir = self.out_dir + "/" + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)+"_"+str(current_time.second)
        self.output_dir = self.output_dir + "_A_" + "{:.2e}".format(self.A) +"_l_" + "{:.2e}".format(self.l)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        # self.record_dir = self.out_dir + "/Record_params"
        # if not os.path.exists(self.record_dir):
        #     os.makedirs(self.record_dir)
        # self.record = self.record_dir + "/record.dat"
#===============================================================================
    def write_sim_params(self):
        with open(self.output_dir+'/sim_params.txt','w') as f:
            attrs= vars(self)
            for param in attrs:
                f.write('{} {}\n'.format(param,attrs[param]))
#===============================================================================
    def make_output_file(self):
        self.output_file= self.output_dir+'/'+'output.out'
        with open(self.output_file, 'w') as f:
            pass
    def copy_anim_script(self):
        subprocess.call('cp {} {}/'.format(self.animscript,self.output_dir), shell=True)
#===============================================================================
    def write_slurm_script(self):
        with open('{}/run.slurm'.format(self.output_dir), 'w') as f:
            f.write("#!/bin/bash -l \n")
            f.write("#SBATCH --job-name=ESGB_run \n")
            f.write("#SBATCH --partition=GravityTheory \n")
            f.write("#SBATCH --mail-type=FAIL \n")
            f.write("#SBATCH --mail-user=ah30@illinois.edu \n")
            f.write("#SBATCH --time=80:00:00 \n")
            f.write("#SBATCH --nodes=1 \n")
            f.write("#SBATCH --ntasks-per-node=1 \n")
            f.write('#SBATCH --mem={}MB\t\t# memory in MB\n'.format(self.memory))
            f.write("#SBATCH --output={}/output.out\n".format(self.output_dir))
            f.write("#SBATCH --error={}/slurm_err.err\n".format(self.output_dir))
            f.write("module load gcc/7.2.0 \n")
            # f.write("conda activate esgb \n")
            f.write("{}/bin/default.run {}".format(self.home_dir,self.output_dir))
#===============================================================================
#===============================================================================
    def write_record(self, line):
        with open(self.record, 'a') as f:
            f.write(line + "\n")
#===============================================================================
    def amplitude_search(self, l, Amp_range, run_type):
        with open(self.record, 'w') as f:
            f.write('run_type = {}\n'.format(run_type))
            f.write('l = {}\n'.format(l))
            f.write('A_low = {}\nA_high = {}\n'.format(Amp_range[0], Amp_range[1]))

        if run_type == "flat_space_to_naked_elliptic":
            A_low = Amp_range[0]
            A_high = Amp_range[1]
            while((A_high - A_low)>1e-3):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if line.startswith("NaN"):
                                self.write_record("NaN; Amp = {}".format(val))
                                A_high = val
                                done = True

                            elif line.startswith("naked"):
                                self.write_record("naked_elliptic_region; Amp = {}".format(val))
                                A_high = val
                                done = True


                            elif line.startswith("exit_code_0"):
                                self.write_record("flat_space; Amp = {}".format(val))
                                A_low = val
                                done = True
                                
        if run_type == "naked_elliptic_to_bh":

            A_low = Amp_range[0]
            A_high = Amp_range[1]
            while((A_high - A_low)>1e-3):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if line.startswith("NaN"):
                                self.write_record("NaN; Amp = {}".format(val))
                                A_low = val
                                done = True

                            elif line.startswith("naked"):
                                self.write_record("naked_elliptic_region; Amp = {}".format(val))
                                A_low = val
                                done = True


                            elif line.startswith("exit_code_1"):
                                self.write_record("bh; Amp = {}".format(val))
                                A_high = val
                                done = True



#===============================================================================
    def launch(self):
        self.make_output_dir()
        if self.slurm == True:
            self.write_sim_params()
            self.copy_anim_script()
            self.write_slurm_script()
            subprocess.call('sbatch {}/run.slurm'.format(self.output_dir),shell = True)
        else:
            self.make_output_file()
            self.write_sim_params()
            self.copy_anim_script()
            subprocess.call('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
        # subprocess.Popen('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
