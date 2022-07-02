#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
import os
#===============================================================================
Amps = np.linspace(1e-4,1,50)
ls = np.linspace(0.01,1,10)
# Amps = [1]
# ls = [1]
input_data = []
for j in range(len(Amps)):
    for l in range(len(ls)):
        input_data.append([ls[l],Amps[j]])


sim = Sim()
# sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory"
sim.animscript = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/Animation-Script.ipynb"
sim.nx = 5000
sim.nt = 5000
sim.save_steps = 100
sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_{}".format(sim.nx,sim.nt)
if not os.path.exists(sim.out_dir):
    os.makedirs(sim.out_dir)
# sim.l = 0.1
sim.initial_mass = 0.
sim.rl = 8.
sim.ru =12.

def launch_sim(vals):
    l_val = vals[0]
    Amp = vals[1]
    sim.l = l_val
    sim.A = Amp
    sim.launch()


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
