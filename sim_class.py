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
    def write_record(self, line):
        with open(self.record, 'a') as f:
            f.write(line + "\n")
#===============================================================================

#===============================================================================
    def compute_initial_mass(self,l,A):
        val = self.collapse_and_bh
        self.l = l
        self.A = A
        self.collapse_and_bh = 0
        self.launch()
        done = False
        while not done:
            time.sleep(2)
            with open(self.output_dir + "/output.out") as f:
                for line in f:
                    if line.startswith("Initial MS_mass = "):
                        index = len("Initial MS_mass = ")
                        M_init = float(line[index:])
                        done = True
        subprocess.call('\nrm -r {}'.format(self.output_dir),shell = True)
        # print("output removed.")
        self.collapse_and_bh = val

        return M_init
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
#===============================================================================
    def amplitude_search(self, l, Amp_range, run_type, tol = 1e-3):
        A_low = Amp_range[0]
        A_high = Amp_range[1]
        M_init = 0
        M_low = self.compute_initial_mass(l,A_low)
        M_high = self.compute_initial_mass(l,A_high)
        with open(self.record, 'w') as f:
            f.write('Bisection search\n')
            f.write('run_type = {}\n'.format(run_type))
            f.write('rl = {}; ru = {}\n'.format(self.rl,self.ru))
            f.write('nx = {}\nnt = {}\n'.format(self.nx,self.nt))
            f.write('l = {}\n'.format(l))
            f.write('A_low = {}; M_init = {}\nA_high = {}; M_init = {}\n'.format(Amp_range[0],M_low, Amp_range[1],M_high))


        if run_type == "flat_space_to_naked_elliptic":

            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                counter = 0
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if(counter==0):
                                if line.startswith("Initial MS_mass = "):
                                    index = len("Initial MS_mass = ")
                                    M_init = float(line[index:])
                                    done = False
                                    counter +=1
                            else:
                                if line.startswith("NaN"):
                                    self.write_record("NaN; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True

                                elif line.startswith("naked"):
                                    self.write_record("naked_elliptic_region; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True


                                elif line.startswith("exit_code_0"):
                                    self.write_record("flat_space; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True

                                elif line.startswith("exit_code_1"):
                                    self.write_record("Problem with run, black hole formed.")
                                    A_high = val
                                    M_high = M_init
                                    done = True


            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; flat_space".format(A_low,M_low))
            self.write_record("A_high = {}; M_init = {}; naked_elliptic_region".format(A_high,M_high))

        elif run_type == "naked_elliptic_to_bh":

            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                counter = 0
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if(counter ==0):
                                if line.startswith("Initial MS_mass = "):
                                    index = len("Initial MS_mass = ")
                                    M_init = float(line[index:])
                                    done = False
                                    counter += 1
                            else:
                                if line.startswith("NaN"):
                                    self.write_record("NaN; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True

                                elif line.startswith("naked"):
                                    self.write_record("naked_elliptic_region; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True


                                elif line.startswith("exit_code_1"):
                                    self.write_record("bh; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True
                                elif line.startswith("exit_code_0"):
                                    self.write_record("Problem with run, flat space formed.")
                                    A_low = val
                                    done = True

            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; naked_elliptic_region".format(A_low, M_low))
            self.write_record("A_high = {}; M_init = {}; bh".format(A_high, M_high))

        elif run_type == "flat_space_fs_to_bh":

            assert l== 0, "l must be 0."
            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                counter = 0
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if counter ==0 :
                                if line.startswith("Initial MS_mass = "):
                                    index = len("Initial MS_mass = ")
                                    M_init = float(line[index:])
                                    done = False
                                    counter +=1
                            else:
                                if line.startswith("NaN"):
                                    self.write_record("NaN; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True

                                elif line.startswith("naked"):
                                    self.write_record("naked_elliptic_region; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True


                                elif line.startswith("exit_code_1"):
                                    self.write_record("bh; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True
                                elif line.startswith("exit_code_0"):
                                    self.write_record("fs; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True

            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; flat_space".format(A_low,M_low))
            self.write_record("A_high = {}; M_init = {}; bh".format(A_high, M_high))

        elif run_type == "collapse_to_bh":

            self.collapse_and_bh = 0
            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.l = l
                self.A = val
                self.launch()
                done = False
                counter = 0
                while not done:
                    time.sleep(30)
                    with open(self.output_dir + "/output.out") as f:
                        for line in f:
                            if counter ==0:
                                if line.startswith("Initial MS_mass = "):
                                    index = len("Initial MS_mass = ")
                                    M_init = float(line[index:])
                                    done = False
                                    counter +=1
                            else:
                                if line.startswith("BH formation at t=0."):
                                    self.write_record("BH_formation; Amp = {}; M_init = {}".format(val,M_init))
                                    A_high = val
                                    M_high = M_init
                                    done = True

                                elif line.startswith("No BH formation at t=0."):
                                    self.write_record("collapse; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True


            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; collapse".format(A_low, M_low))
            self.write_record("A_high = {}; M_init = {}; bh".format(A_high, M_high))
