#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
# Amps = np.linspace(1e-2,1,100)
# ls = np.linspace(1e-1,10,100)

Amps = np.array([0.3,5])
ls = np.array([1])

input_data = []

for j in range(len(Amps)):
    for l in range(len(ls)):
        input_data.append([ls[l],Amps[j]])

current_time = datetime.now()
sim = Sim()
sim.slurm = False
# sim.animscript = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/Animation-Script.ipynb"
sim.animscript = "./Animation-Script.ipynb"
sim.nx = 20000
sim.nt = 20000
sim.save_steps = int(sim.nt/10)
# sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_".format(sim.nx,sim.nt)+ current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
sim.out_dir = "./output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_{}_".format(sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)

if not os.path.exists(sim.out_dir):
    os.makedirs(sim.out_dir)

sim.initial_mass = 0
sim.exc_i = 0
sim.rl = 8.
sim.ru =12.
#====================================================
run_params = sim.out_dir + "/Run_params"
if not os.path.exists(run_params):
    os.makedirs(run_params)

with open(run_params + "/run_params.dat", "w" ) as f:
    f.write("Total number of runs = {}\n".format(len(input_data)))
    f.write("nx = {} \n".format(sim.nx))
    f.write("nt = {} \n".format(sim.nt))
    f.write("save_steps = {} \n".format(sim.save_steps))


np.savetxt(run_params + "/ls.dat" , ls)
np.savetxt(run_params + "/Amps.dat", Amps)
#===================================================
def launch_sim(vals):
    l_val = vals[0]
    Amp = vals[1]
    sim.l = l_val
    sim.A = Amp
    sim.launch()

if sim.slurm == True:
    sim.memory = '30'
    # for j in range(len(Amps)):
    #     # launch_sim([ls[j],Amps[j]])
    #     for l in range(len(ls)):
    #         launch_sim([ls[l],Amps[j]])
    for j in range(len(input_data)):
        launch_sim(input_data[j])
else:
    if __name__ == '__main__':
        pool_nums = 7
        t_start = time.time()
        print("Starting multiprocessing pool..")
        pool = Pool(pool_nums)
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
        print("Data saved at:{}".format(sim.out_dir))

        with open(run_params + "/run_params.dat", "w" ) as f:
            f.write("Number of pools = {}\n".format(pool_nums))
            f.write("Finished process. \nTime = {} s\n".format(t_end- t_start))
            f.write("Data saved at:{} \n".format(sim.out_dir))

        # print("Ending caffeinate...")
        # subprocess.call("killall caffeinate" , shell = True)
        # print("Done. Exit.")
#===============================================================================
