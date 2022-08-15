#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
#theory = "shift_symm"
#theory = "gaussian"
home_path = "."
#home_path = "/home/ah30/scratch/code-f-phi"
Amps = np.array([0.12])
ls = np.array([1.])
lexp = np.array([0.])
mu = np.array([0.])

out_path = home_path+ "/output/Phase-Space/Runs_all"
#===============================================================================
input_data = []

# for j in range(len(Amps)):
#     for l in range(len(ls)):
#         input_data.append([ls[l],Amps[j]])
# input_data = [[0.0184375,0.,1.,12.],
# [0.0203125,0.,0.9,12.],
# [0.02171875,0.,0.8,12.],
# [0.02328125,0.,0.7,12.],
# [0.025625,0.,0.6,12.],
# [0.027,0.,0.5,12.],
# [0.030,0.,0.4,12.],
# [0.033,0.,0.3,12.]]
input_data = [
[0.0367187500000000,0.1,0.,0.]]
input_data = np.array(input_data)
current_time = datetime.now()
sim = Sim()
sim.slurm = False
sim.cluster = False
sim.write_runs = True
sim.animscript = home_path +"/Animation-Script.ipynb"
sim.cl = 100.0
sim.nx = 5000
sim.nt = 5000
sim.save_steps = int(sim.nt/10)
sim.initial_mass = 0
if(sim.initial_mass == 0):
    sim.exc_i = 0
else:
    sim.exc_i = 3
sim.exc_i = 0
sim.rl = 8.
sim.ru =12.
sim.collapse_and_bh = 1;
sim.search =False
#===============================================================================
if sim.search == True:
    sim.out_dir = out_path+"/Search/Search_rl_{}_ru_{}/Run_nx_{}_nt_{}_".format(sim.rl,sim.ru,sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
else:
    sim.out_dir = out_path + "/Runs/Runs_rl_{}_ru_{}/Run_nx_{}_nt_{}_".format(sim.rl,sim.ru,sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)

if not os.path.exists(sim.out_dir):
    os.makedirs(sim.out_dir)

#===============================================================================
run_params = sim.out_dir + "/Run_params"
if not os.path.exists(run_params):
    os.makedirs(run_params)

if sim.search == False:
    with open(run_params + "/run_params.dat", "w" ) as f:
        f.write("Total number of runs = {}\n".format(len(input_data)))
        f.write("nx = {} \n".format(sim.nx))
        f.write("nt = {} \n".format(sim.nt))
        f.write("save_steps = {} \n".format(sim.save_steps))

    np.savetxt(run_params + "/Amps.dat", input_data[:,0])
    np.savetxt(run_params + "/ls.dat" , input_data[:,1])
    np.savetxt(run_params + "/lexp.dat" , input_data[:,2])
    np.savetxt(run_params + "/mu.dat" , input_data[:,3])
else:
    with open(run_params + "/run_params.dat", "w" ) as f:
        f.write("nx = {} \n".format(sim.nx))
        f.write("nt = {} \n".format(sim.nt))
        f.write("save_steps = {} \n".format(sim.save_steps))

#===============================================================================
def launch_sim(vals):
    Amp = vals[0]
    l_s = vals[1]
    l_exp = vals[2]
    mu_s = vals[3]
    sim.ls = l_s
    sim.lexp = l_exp
    sim.mu = mu_s
    sim.A = Amp
    sim.launch()
#===============================================================================
if sim.search == True:
    data_search = [[0.1,0.2,1e-3,0,0.3,3],
    [0.1,0.2,1e-3,0,0.4,3],
    [0.1,0.2,1e-3,0,0.5,3],
    [0.1,0.2,1e-3,0,0.6,3],
    [0.1,0.2,1e-3,0,0.7,3],
    [0.1,0.2,1e-3,0,0.8,3],
    [0.1,0.2,1e-3,0,0.9,3],
    [0.1,0.2,1e-3,0,1.,3]
      ]



    #["flat_space_to_naked_elliptic","naked_elliptic_to_blackhole","flat_space_fs_to_blackhole","collapse_to_blackhole"]

    run_type = "naked_elliptic_to_blackhole"

    def launch_search(arr):
        Amp_range = [arr[0],arr[1]]
        tol = arr[2]
        ls = arr[3]
        lexp = arr[4]
        mu = arr[5]
        sim.record = run_params + "/record_ls_{}_lexp_{}_mu_{}.dat".format(ls,lexp,mu)
        sim.amplitude_search(ls=ls,lexp = lexp,mu=mu, Amp_range = Amp_range , run_type = run_type, tol = tol)

    #--------------------------------------------------------------------------
    if __name__ == '__main__':
        # print("theory = ",theory)
        if sim.cluster:
            pool_nums = len(data_search)
        else:
            if len(data_search) >=6:
                pool_nums = 6
            else :
                pool_nums = len(data_search)



        print("pool_nums = ", pool_nums)

        t_start = time.time()

        print("Searching amplitude...")
        print("run_type = {}".format(run_type))
        print("nx = ",sim.nx)
        print("nt = ", sim.nt)
        print("save_steps = ", sim.save_steps)
        print(data_search)

        print("Data saved at:{}".format(sim.out_dir))

        print("Starting multiprocessing pool..")

        pool = Pool(pool_nums)

        result = pool.map_async(launch_search, data_search)

        pool.close()
        while True:
            if not result.ready():
                print('We\'re not done yet, %s tasks to go!' % result._number_left)
                time.sleep(60)
            else:
                break

        pool.join()

        t_end = time.time()

        print("Finished process. \nTime = ",t_end-t_start," s")

        with open(run_params + "/run_params.dat", "a" ) as f:
            f.write("Finished process. \nTime = {} s\n".format(t_end- t_start))
            f.write("Data saved at:{} \n".format(sim.out_dir))
else:
    if __name__ == '__main__':
        if len(input_data) >=6:
            pool_nums = 6
        else :
            pool_nums = len(input_data)


        print("pool_nums = ", pool_nums)

        t_start = time.time()
        print("Starting multiprocessing pool..")
        print("Data saved at:{}".format(sim.out_dir))
        print("nx = ",sim.nx)
        print("nt = ", sim.nt)
        print("save_steps = ", sim.save_steps)
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

        with open(run_params + "/run_params.dat", "a" ) as f:
            f.write("Number of pools = {}\n".format(pool_nums))
            f.write("Finished process. \nTime = {} s\n".format(t_end- t_start))
            f.write("Data saved at:{} \n".format(sim.out_dir))
#===============================================================================
