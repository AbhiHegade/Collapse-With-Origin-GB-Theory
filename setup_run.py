#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
theory = "shift_symm"
# theory = "gaussian"
print("theory = ",theory)
Amps = np.array([1e-3])
ls = np.array([2])

if theory == "shift_symm":
    out_path = "./output/Phase-Space/Shift-Symmetric-Theory"
else:
    out_path = "./output/Phase-Space/Gaussian"
#===============================================================================
input_data = []

for j in range(len(Amps)):
    for l in range(len(ls)):
        input_data.append([ls[l],Amps[j]])

current_time = datetime.now()
sim = Sim()
sim.slurm = False
sim.animscript = "./Animation-Script.ipynb"
sim.nx = 10000
sim.nt = 10000
sim.save_steps = int(sim.nt/10)
sim.initial_mass = 0
sim.exc_i = 0
sim.rl = 20.
sim.ru =24.
sim.collapse_and_bh = 1;
sim.search =True
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
        # data_search = []
        # ls = np.linspace(1,2,11)
        #
        # for j in range(len(ls)):
        #     data_search.append([ls[j],1e-3 ,2e-2])
        #
        # data_search = np.array(data_search)
        data_search = np.array([
       [1.000000e-01, 6.562500e-03, 1.050000e-01],
       [1.500000e-01, 4.609375e-03, 7.375000e-02],
       [2.000000e-01, 3.281250e-03, 5.250000e-02],
       [2.500000e-01, 2.375000e-03, 3.800000e-02],
       [3.000000e-01, 1.718750e-03, 2.750000e-02],
       [3.500000e-01, 1.469580e-03, 2.351328e-02],
       [4.000000e-01, 1.000000e-03, 1.600000e-02],
       [4.500000e-01, 9.706550e-04, 1.553048e-02],
       [5.000000e-01, 8.066400e-04, 1.290624e-02],
       [5.500000e-01, 6.860475e-04, 1.097676e-02],
       [6.000000e-01, 5.917975e-04, 9.468760e-03],
       [6.500000e-01, 5.361425e-04, 8.578280e-03],
       [7.000000e-01, 4.375000e-04, 7.000000e-03],
       [7.500000e-01, 3.880325e-04, 6.208520e-03],
       [8.000000e-01, 3.625000e-04, 5.800000e-03],
       [8.500000e-01, 3.041425e-04, 4.866280e-03],
       [9.000000e-01, 2.750000e-04, 4.400000e-03],
       [9.500000e-01, 2.445525e-04, 3.912840e-03],
       [1.000000e+00, 2.257825e-04, 3.612520e-03]])

        tol = 1e-4



        run_type = "flat_space_to_naked_elliptic"

        def launch_search(arr):
            l = arr[0]
            Amp_range = [arr[1],arr[2]]
            sim.record = run_params + "/record_{}.dat".format(l)
            sim.amplitude_search(l=l, Amp_range = Amp_range , run_type = run_type, tol = tol)

        #--------------------------------------------------------------------------
        if __name__ == '__main__':
            if len(data_search) >=6:
                pool_nums = 6
            else :
                pool_nums = len(data_search)



            print("pool_nums = ", pool_nums)

            t_start = time.time()

            print("Searching amplitude...")
            print("run_type = {}".format(run_type))

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

            t_end = time.time()

            print("Finished process. \nTime = ",t_end-t_start," s")

            with open(run_params + "/run_params.dat", "w" ) as f:
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

            with open(run_params + "/run_params.dat", "w" ) as f:
                f.write("Number of pools = {}\n".format(pool_nums))
                f.write("Finished process. \nTime = {} s\n".format(t_end- t_start))
                f.write("Data saved at:{} \n".format(sim.out_dir))

        # print("Ending caffeinate...")
        # subprocess.call("killall caffeinate" , shell = True)
        # print("Done. Exit.")
#===============================================================================
