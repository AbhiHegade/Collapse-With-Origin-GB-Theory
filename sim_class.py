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
        self.output_dir = self.output_dir+"_M_{:.2e}".format(self.initial_mass) + "_A_" + "{:.2e}".format(self.A) +"_ls_" + "{:.2e}".format(self.ls) + "_lexp_" + "{:.2e}".format(self.lexp) + "_mu_" + "{:.2e}".format(self.mu)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
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
    def compute_initial_mass(self,A,ls,lexp,mu):
        val = self.collapse_and_bh
        self.A = A
        self.ls = ls
        self.lexp = lexp
        self.mu = mu
        self.collapse_and_bh = 0
        self.launch()
        time.sleep(2)
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
            if self.write_runs:
                with open("{}/output.out".format(self.output_dir),'r') as f:
                    for line in f:
                        if line.startswith("Initial MS_mass = "):
                            mass_init = float(line[len("Initial MS_mass = "):])
                    run_status = line
                with open("{}/Run_params/run_params.dat".format(self.out_dir), 'a') as f:
                    f.write("A = {}; ls = {}; lexp = {}; mu = {}; M_init = {}; status ={}".format(self.A,self.ls,self.lexp,self.mu,mass_init,run_status))


        # subprocess.Popen('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
#===============================================================================
    def amplitude_search(self, ls,lexp,mu,Amp_range, run_type, tol = 1e-3):
        A_low = Amp_range[0]
        A_high = Amp_range[1]
        M_init = 0
        M_low = self.compute_initial_mass(A_low,ls,lexp,mu)
        M_high = self.compute_initial_mass(A_high,ls,lexp,mu)
        with open(self.record, 'w') as f:
            f.write('Bisection search\n')
            f.write('tol = {}\n'.format(tol))
            f.write('run_type = {}\n'.format(run_type))
            f.write('rl = {}; ru = {}\n'.format(self.rl,self.ru))
            f.write('nx = {}\nnt = {}\n'.format(self.nx,self.nt))
            f.write('ls = {}\n'.format(ls))
            f.write('lexp = {}\n'.format(lexp))
            f.write('mu = {}\n'.format(mu))
            f.write('A_low = {}; M_init = {}\nA_high = {}; M_init = {}\n'.format(Amp_range[0],M_low, Amp_range[1],M_high))

        if run_type == "flat_space_to_naked_elliptic":

            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.ls = ls
                self.lexp = lexp
                self.mu = mu
                self.A = val
                self.launch()
                # time.sleep(60)
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
                                    self.write_record("NaN; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_high = val
                                    M_high = M_init
                                    done = True

                                elif line.startswith("naked"):
                                    self.write_record("naked_elliptic_region; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_high = val
                                    M_high = M_init
                                    done = True


                                elif line.startswith("exit_code_0"):
                                    self.write_record("flat_space; Amp = {}; M_init = {}".format(val,M_init))
                                    A_low = val
                                    M_low = M_init
                                    done = True

                                elif line.startswith("exit_code_1"):
                                    self.write_record("Problem with run, black hole formed; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_high = val
                                    M_high = M_init
                                    done = True


            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; flat_space".format(A_low,M_low))
            self.write_record("A_high = {}; M_init = {}; naked_elliptic_region".format(A_high,M_high))

        elif run_type == "naked_elliptic_to_blackhole":

            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.ls = ls
                self.lexp = lexp
                self.mu = mu
                self.A = val
                self.launch()
                #time.sleep(60)
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
                                    self.write_record("NaN; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_low = val
                                    M_low = M_init
                                    done = True

                                elif line.startswith("naked"):
                                    self.write_record("naked_elliptic_region; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_low = val
                                    M_low = M_init
                                    done = True


                                elif line.startswith("exit_code_1"):
                                    self.write_record("bh; Amp = {}; M_init = {}; status = {}".format(val,M_init,line))
                                    A_high = val
                                    M_high = M_init
                                    done = True
                                elif line.startswith("exit_code_0"):
                                    self.write_record("Problem with run, flat space formed; status = {}".format(line))
                                    A_low = val
                                    done = True

            self.write_record("Run finished.")
            self.write_record("A_low = {}; M_init = {}; naked_elliptic_region".format(A_low, M_low))
            self.write_record("A_high = {}; M_init = {}; bh".format(A_high, M_high))

        elif run_type == "flat_space_fs_to_blackhole":

            assert ls== 0, "ls must be 0."
            assert lexp== 0, "lexp must be 0."
            assert mu== 0, "mu must be 0."
            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.ls = ls
                self.lexp = lexp
                self.mu = mu
                self.A = val
                self.launch()
                time.sleep(60)
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

        elif run_type == "collapse_to_blackhole":

            self.collapse_and_bh = 0
            while((A_high - A_low)>tol):
                val = (A_high + A_low)/2
                self.ls = ls
                self.lexp = lexp
                self.mu = mu
                self.A = val
                self.launch()
                # time.sleep(60)
                done = False
                counter = 0
                while not done:
                    time.sleep(3)
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
#===============================================================================
    def mass_search(self, ls, lexp, mu, mass_range, tol = 1e-3):
        M_low = mass_range[0]
        M_high = mass_range[1]
        run_type = "black_hole_mass_search"
        with open(self.record, 'w') as f:
            f.write('Bisection search\n')
            f.write('run_type = {}\n'.format(run_type))
            f.write('tol = {}\n'.format(tol))
            f.write('A = {}; rl = {}; ru = {}\n'.format(self.A,self.rl,self.ru))
            f.write('nx = {}\nnt = {}\n'.format(self.nx,self.nt))
            f.write('ls = {}\n'.format(ls))
            f.write('lexp = {}\n'.format(lexp))
            f.write('mu = {}\n'.format(mu))
            f.write('tol = {}\n'.format(tol))
            f.write('M_low = {}; \nM_high = {};\n'.format(M_low,M_high))
        while((M_high - M_low)>tol):
            val = (M_high + M_low)/2
            self.ls = ls
            self.lexp = lexp
            self.mu = mu
            self.initial_mass = val
            dx = self.cl/self.nx
            ratio = 0.5
            self.exc_i = int((ratio/dx)*((2*self.cl*self.initial_mass)/(self.cl + 2*self.initial_mass)))
            self.launch()
            time.sleep(2)
            done = False
            while not done:
                time.sleep(30)
                with open(self.output_dir + "/output.out") as f:
                    for line in f:
                        if line.startswith("NaN"):
                            self.write_record("NaN; Mass = {}; status = {}".format(val,line))
                            M_low = val
                            done = True

                        elif line.startswith("naked"):
                            self.write_record("naked_elliptic_region; Mass = {}; status = {}".format(val,line))
                            M_low = val
                            done = True

                        elif line.startswith("exit_code_1"):
                            self.write_record("bh; Mass = {}; status = {}".format(val,line))
                            M_high = val
                            done = True

                        elif line.startswith("exit_code_0"):
                            self.write_record("Problem with run with run flat space formed; Mass = {}; status = {}".format(val,line))
                            M_low = val
                            done = True



        self.write_record("Run finished.")
        self.write_record("M_low = {}; naked_elliptic_region".format(M_low))
        self.write_record("M_high = {}; bh".format(M_high))
