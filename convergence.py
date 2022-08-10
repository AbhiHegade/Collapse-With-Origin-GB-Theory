#!/usr/bin/env python
# coding: utf-8
import os, sys, time, shutil, subprocess
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os

class Convg:
    def __init__(self):
        self.home_dir = str(os.getcwd())

    def make_output_dir(self,level):
        current_time = datetime.now()
        self.output_dir = self.out_dir + "/" + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)+"_"+str(current_time.second)
        self.output_dir = self.output_dir + "_M_{:.2e}".format(self.initial_mass) + "_A_" + "{:.2e}".format(self.A) +"_ls_" + "{:.2e}".format(self.ls) + "_lexp_" + "{:.2e}".format(self.lexp) + "_mu_" + "{:.2e}".format(self.mu)
        self.output_dir = self.output_dir  + "/{}".format(level)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if level == "h":
            dir_copy_to = "/".join(self.output_dir.split("/")[:-1])
            subprocess.call('cp {} {}/'.format(self.animscript,dir_copy_to), shell=True)

    def make_output_file(self):
        self.output_file= self.output_dir+'/'+'output.out'
        with open(self.output_file, 'w') as f:
            pass

    def write_sim_params(self):
        with open(self.output_dir+'/sim_params.txt','w') as f:
            attrs= vars(self)
            for param in attrs:
                f.write('{} {}\n'.format(param,attrs[param]))

    def copy_anim_script(self):
        subprocess.call('cp {} {}/'.format(self.animscript,self.output_dir), shell=True)

    def launch(self,level):
        self.make_output_dir(level)
        self.make_output_file()
        self.copy_anim_script()
        self.write_sim_params()
        subprocess.call('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)

def generate_convg_arr(arr):
    rnums = arr.shape[0]
    return_arr = np.empty(shape = (3*rnums,5))
    for j in range(len(arr)):
        Amp,ls,lexp,mu = arr[j]
        j1 = 3*j
        return_arr[j1] = np.array([Amp,ls,lexp,mu,4])
        return_arr[j1+1] = np.array([Amp,ls,lexp,mu,2])
        return_arr[j1+2] = np.array([Amp,ls,lexp,mu,1])
    return return_arr
def generate_convg_arr_mass(arr):
    rnums = arr.shape[0]
    return_arr = np.empty(shape = (3*rnums,5))
    for j in range(len(arr)):
        mass,ls,lexp,mu = arr[j]
        j1 = 3*j
        return_arr[j1] = np.array([mass,ls,lexp,mu,4])
        return_arr[j1+1] = np.array([mass,ls,lexp,mu,2])
        return_arr[j1+2] = np.array([mass,ls,lexp,mu,1])
    return return_arr

def get_arr(arr,r_type):
    arr = np.array(arr)
    if r_type == "normal":
        arr = generate_convg_arr(arr)
    elif r_type == "mass":
        arr = generate_convg_arr_mass(arr)

    return arr
#===============================================================================
home_path = "."
#home_path = "/home/ah30/scratch/code-f-phi"
Amps = np.array([0.12])
ls = np.array([0.5])
lexp = np.array([0.])
mu = np.array([0.])

out_path = home_path+ "/output/Phase-Space/Convg"
#===============================================================================
input_data = []
# input_data = [[0.16,0.5,0,0],
# [0.02,0.5,0,0]]
# input_data = [
# [0.01,0.5,0,0]]
# input_data_mass = [[2.,0.5,0,0], [0.45,0.,0.5,3]]
# input_data_mass = [[3.,0.1,0,0],[3.,0.5,0,0]]
# input_data_mass = get_arr(input_data_mass,"mass")
input_data = [[0.16,0.5,0,0]]
input_data = get_arr(input_data,"normal")
current_time = datetime.now()
sim = Convg()
sim.animscript = home_path +"/Animation-Script.ipynb"
sim.mass_run = False
sim.cl = 100.0
sim.initial_mass = 0
if(sim.initial_mass == 0):
    sim.exc_i = 0
else:
    sim.exc_i = 3
sim.animscript = home_path +"/Convergence-Analysis.ipynb"
nx = 4000
nt = 4000
sim.exc_i = 0
sim.rl = 8.
sim.ru =12.
sim.collapse_and_bh = 1;
#===============================================================================
sim.out_dir = out_path+"/Convg_rl_{}_ru_{}/Run_nx_{}_nt_{}_".format(sim.rl,sim.ru,nx,nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
if not os.path.exists(sim.out_dir):
    os.makedirs(sim.out_dir)

#===============================================================================
def launch_sim(vals):
    Amp = vals[0]
    l_s = vals[1]
    l_exp = vals[2]
    mu_s = vals[3]
    level = vals[4]
    sim.ls = l_s
    sim.lexp = l_exp
    sim.mu = mu_s
    sim.A = Amp
    if level == 4:
        sim.nx = nx
        sim.nt = nt
        sim.save_steps = int(sim.nt/100)
        level_str = "4h"
    elif level == 2:
        sim.nx = 2*nx
        sim.nt = 2*nt
        sim.save_steps = int(sim.nt/100)
        level_str = "2h"
    elif level ==1:
        sim.nx = 4*nx
        sim.nt = 4*nt
        sim.save_steps = int(sim.nt/100)
        level_str = "h"
    sim.launch(level_str)
#===============================================================================
def launch_sim_mass(vals):
    mass = vals[0]
    l_s = vals[1]
    l_exp = vals[2]
    mu_s = vals[3]
    level = vals[4]
    sim.ls = l_s
    sim.lexp = l_exp
    sim.mu = mu_s
    sim.A = 1e-3
    if level == 4:
        sim.nx = nx
        sim.nt = nt
        sim.save_steps = int(sim.nt/100)
        sim.initial_mass = mass
        dx = sim.cl/sim.nx
        ratio = 0.5
        sim.exc_i = int((ratio/dx)*((4*sim.cl*sim.initial_mass)/(sim.cl + 4*sim.initial_mass)))
        level_str = "4h"
    elif level == 2:
        sim.nx = 2*nx
        sim.nt = 2*nt
        sim.save_steps = int(sim.nt/100)
        sim.initial_mass = mass
        dx = sim.cl/sim.nx
        ratio = 0.5
        sim.exc_i = int((ratio/dx)*((4*sim.cl*sim.initial_mass)/(sim.cl + 4*sim.initial_mass)))
        level_str = "2h"
    elif level ==1:
        sim.nx = 4*nx
        sim.nt = 4*nt
        sim.save_steps = int(sim.nt/100)
        sim.initial_mass = mass
        dx = sim.cl/sim.nx
        ratio = 0.5
        sim.exc_i = int((ratio/dx)*((4*sim.cl*sim.initial_mass)/(sim.cl + 4*sim.initial_mass)))
        level_str = "h"
    sim.launch(level_str)
#===============================================================================
if __name__ == '__main__':
    if sim.mass_run == True:
        if len(input_data_mass) >=6:
            pool_nums = 6
        else :
            pool_nums = len(input_data_mass)


        print("pool_nums = ", pool_nums)

        t_start = time.time()
        print("Starting multiprocessing pool..")
        print("Convergence Run")
        print("Data saved at:{}".format(sim.out_dir))
        print("nx = ",nx)
        print("nt = ", nt)
        pool = Pool(pool_nums)
        result = pool.map_async(launch_sim_mass, input_data_mass)
        pool.close()

        while True:
            if not result.ready():
                print('We\'re not done yet, %s tasks to go!' % result._number_left)
                time.sleep(20)
            else:
                break

        pool.join()

        t_end = time.time()
        print("Finished process. \nTime = ",t_end-t_start," s")
        print("Data saved at:{}".format(sim.out_dir))

    else :
        if len(input_data) >=6:
            pool_nums = 6
        else :
            pool_nums = len(input_data)


        print("pool_nums = ", pool_nums)

        t_start = time.time()
        print("Starting multiprocessing pool..")
        print("Convergence Run")
        print("Data saved at:{}".format(sim.out_dir))
        print("nx = ",nx)
        print("nt = ", nt)
        pool = Pool(pool_nums)
        result = pool.map_async(launch_sim, input_data)
        pool.close()

        while True:
            if not result.ready():
                print('We\'re not done yet, %s tasks to go!' % result._number_left)
                time.sleep(20)
            else:
                break

        pool.join()

        t_end = time.time()
        print("Finished process. \nTime = ",t_end-t_start," s")
        print("Data saved at:{}".format(sim.out_dir))
