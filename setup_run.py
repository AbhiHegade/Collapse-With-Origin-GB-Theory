#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
theory = "gaussian"
Amps = np.array([0.24])
ls = np.array([0.1])

if theory == "shift_symm":
    out_path = "./output/Phase-Space/Shift-Symmetric-Theory"
else:
    out_path = "./output/Phase-Space/Gaussian"

input_data = []

for j in range(len(Amps)):
    for l in range(len(ls)):
        input_data.append([ls[l],Amps[j]])

# input_data = np.array(input_data)[400:1000]
current_time = datetime.now()
sim = Sim()
sim.slurm = False
# sim.animscript = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/Animation-Script.ipynb"
sim.animscript = "./Animation-Script.ipynb"
sim.nx = 5000
sim.nt = 5000
sim.save_steps = int(sim.nt/10)
sim.search =False
# sim.out_dir = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Phase-Space/Shift-Symmetric-Theory/Run_nx_{}_nt_".format(sim.nx,sim.nt)+ current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
if sim.search == True:
    sim.out_dir = out_path+"/Search/Run_nx_{}_nt_{}_".format(sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
else:
    sim.out_dir = out_path + "/Runs/Run_nx_{}_nt_{}_".format(sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)

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
if sim.search == False:
    with open(run_params + "/run_params.dat", "w" ) as f:
        f.write("Total number of runs = {}\n".format(len(input_data)))
        f.write("nx = {} \n".format(sim.nx))
        f.write("nt = {} \n".format(sim.nt))
        f.write("save_steps = {} \n".format(sim.save_steps))


    np.savetxt(run_params + "/ls.dat" , ls)
    np.savetxt(run_params + "/Amps.dat", Amps)
else:
#===================================================
    with open(run_params + "/run_params.dat", "w" ) as f:
        f.write("nx = {} \n".format(sim.nx))
        f.write("nt = {} \n".format(sim.nt))
        f.write("save_steps = {} \n".format(sim.save_steps))

#===================================================
def launch_sim(vals):
    l_val = vals[0]
    Amp = vals[1]
    sim.l = l_val
    sim.A = Amp
    sim.launch()

if sim.slurm == True:
    sim.memory = '30'
    if sim.search == True:
        l = 0.1
        Amp_range = [0.025,0.0275]
        run_type = "flat_space_to_naked_elliptic"
        sim.amplitude_search(l=l, Amp_range = Amp_range , run_type = run_type)
    else:
        for j in range(len(input_data)):
            launch_sim(input_data[j])

else:
    if sim.search == True:
        # data_search = np.array([[0.5,1e-3,6e-3], [0.6,1e-3,6e-3],
        # [0.7, 1e-3, 5e-3], [0.8, 9e-4, 2e-3], [0.9, 8e-4,2e-3],
        # [1,7e-4,2e-3]])
        data_search = np.array([[0.1, 0.2, 0.4], [0.15,0.2,0.4], [0.2,0.2,0.4],[0.25,0.2,0.4], [0.3,0.2,0.4 ]])
        # data_search = np.array([[0.1, 0.25, 0.35]])
        # data_search = []
        # for x in np.linspace(0.1,1 ,19)[5:-2]:
        #     data_search.append([x,0.45, 0.7])
        #
        # data_search = np.array(data_search)
        tol = 1e-3
        pool_nums = 5
        run_type = "naked_elliptic_to_bh"
        def launch_search(arr):
            l = arr[0]
            Amp_range = [arr[1],arr[2]]
            sim.record = run_params + "/record_{}.dat".format(l)
            sim.amplitude_search(l=l, Amp_range = Amp_range , run_type = run_type, tol = tol)

        #--------------------------------------------------------------------------
        if __name__ == '__main__':

            t_start = time.time()
            print("Searching amplitude...")
            print("run_type = {}".format(run_type))
            # print("l = {}".format(data_search[:,0]))
            # print("Amp_range = {}".format(Amp_range))
            print(data_search)
            print("Data saved at:{}".format(sim.out_dir))

            print("Starting multiprocessing pool..")
            pool = Pool(pool_nums)
            result = pool.map_async(launch_search, data_search)

            pool.close()
            while True:
                if not result.ready():
                    print('We\'re not done yet, %s tasks to go!' % result._number_left)
                    time.sleep(20)
                else:
                    break

            pool.join()
            # for xx in data_search:
            #     launch_search(xx)
            t_end = time.time()

            print("Finished process. \nTime = ",t_end-t_start," s")

            with open(run_params + "/run_params.dat", "w" ) as f:
                f.write("Finished process. \nTime = {} s\n".format(t_end- t_start))
                f.write("Data saved at:{} \n".format(sim.out_dir))
    else:
        if __name__ == '__main__':
            pool_nums = 1
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
