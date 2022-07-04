#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
Amps = np.array([0.4694408163265306,1.0,0.9591877551020408,0.5918775510204082,0.9387816326530612,0.9795938775510205,0.897969387755102,0.9183755102040816])
# ls = np.linspace(0.01,1,1)
ls = np.array([1])
# Amps = [1]
# ls = [1]
input_data = []
for j in range(len(Amps)):
    for l in range(len(ls)):
        input_data.append([ls[l],Amps[j]])

current_time = datetime.now()
sim = Sim()
sim.slurm = False
# sim.animscript = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/Animation-Script.ipynb"
sim.animscript = "./Animation-Script.ipynb"
sim.nx = 5000
sim.nt = 5000
sim.save_steps = 100
# sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_".format(sim.nx,sim.nt)+ current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
sim.out_dir = "./output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_{}_".format(sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
if not os.path.exists(sim.out_dir):
    os.makedirs(sim.out_dir)
sim.initial_mass = 0.
sim.rl = 8.
sim.ru =12.

def launch_sim(vals):
    l_val = vals[0]
    Amp = vals[1]
    sim.l = l_val
    sim.A = Amp
    sim.launch()

if sim.slurm == True:
    sim.memory = '30'
    for j in range(len(Amps)):
        # launch_sim([ls[j],Amps[j]])
        for l in range(len(ls)):
            launch_sim([ls[l],Amps[j]])
else:
    if __name__ == '__main__':
        t_start = time.time()
        print("Starting multiprocessing pool..")
        pool = Pool(7)
        result = pool.map_async(launch_sim, input_data)
        pool.close()

        while True:
            if not result.ready():
                print('We\'re not done yet, %s tasks to go!' % result._number_left)
                time.sleep(20)
            else:
                break
        # pool.close()
        # pool.join()
        pool.join()

        t_end = time.time()
        print("Finished process. \nTime = ",t_end-t_start," s")
#===============================================================================
