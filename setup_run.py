#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
#===============================================================================
Amps = np.linspace(1e-4,1,10)
ls = np.linspace(1e-4,1,10)
input_data = []
for l in range(len(ls)):
    for j in range(len(Amps)):
        input_data.append([Amps[j],ls[l]])

sim = Sim()
sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory"
sim.nx = 4000
sim.nt = 4000
sim.save_steps = 100
sim.l = 0.1
sim.initial_mass = 0.
sim.rl = 8.
sim.ru =12.

def launch_sim(vals):
    Amp = vals[0]
    l_val = vals[1]
    sim.A = Amp
    sim.l = l_val
    sim.launch()


if __name__ == '__main__':
    t_start = time.time()
    print("Starting multiprocessing pool..")
    pool = Pool(5)
    pool.map(launch_sim, input_data)
    pool.close()
    pool.join()


    t_end = time.time()
    print("Finished process. \nTime = ",t_end-t_start," s")
#===============================================================================
