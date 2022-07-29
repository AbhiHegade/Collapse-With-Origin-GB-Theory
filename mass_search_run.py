#!/usr/bin/env python
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
import os
#===============================================================================
#theory = "shift_symm"
theory = "gaussian"
#home_path = "."
home_path = "/home/ah30/scratch/code-f-phi"
Ms = np.array([0.9,0.8,1.])
ls = np.array([1.])

if theory == "shift_symm":
    out_path = home_path+ "/output/Phase-Space/Shift-Symmetric-Theory"
else:
    out_path = home_path+ "/output/Phase-Space/Gaussian"
#===============================================================================

#input_data = []

# for j in range(len(Ms)):
#     for l in range(len(ls)):
#         assert (ls[l]>0), "l must be greater than zero."
#         input_data.append([ls[l],Ms[j]])

input_data  = [[1,0.98],[1,1.15]]
input_data = np.array(input_data)
current_time = datetime.now()
sim = Sim()
sim.slurm = False
sim.write_runs = False
sim.animscript = home_path+ "/Animation-Script.ipynb"
sim.cl = 100.0
sim.nx = 4000
sim.nt = 160000
sim.save_steps = int(sim.nt/10)
sim.initial_mass = 1
if(sim.initial_mass == 0):
    sim.exc_i = 0
else:
    sim.exc_i = 3
sim.A = 1e-2
sim.rl = 8.
sim.ru =12.
sim.collapse_and_bh = 1;
sim.search =True
#===============================================================================
if sim.search == True:
    sim.out_dir = out_path+"/Search/Search_Mass_rl_{}_ru_{}/Run_nx_{}_nt_{}_".format(sim.rl,sim.ru,sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)
else:
    sim.out_dir = out_path + "/Runs/Runs_Mass_rl_{}_ru_{}/Run_nx_{}_nt_{}_".format(sim.rl,sim.ru,sim.nx,sim.nt) + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_"+str(current_time.minute)

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


    np.savetxt(run_params + "/ls.dat" , input_data[:,0])
    np.savetxt(run_params + "/masses.dat", input_data[:,1])
else:
#===================================================
    with open(run_params + "/run_params.dat", "w" ) as f:
        f.write("nx = {} \n".format(sim.nx))
        f.write("nt = {} \n".format(sim.nt))
        f.write("save_steps = {} \n".format(sim.save_steps))

#===================================================
def launch_sim(vals):
    l_val = vals[0]
    Mass_val = vals[1]
    sim.l = l_val
    sim.initial_mass = Mass_val
    dx = sim.cl/sim.nx
    ratio = 0.5
    sim.exc_i = int((ratio/dx)*((4*sim.cl*sim.initial_mass)/(sim.cl + 4*sim.initial_mass)))
    sim.launch()

if sim.slurm == True:
    sim.memory = '30'
    if sim.search == True:
        l = 0.5
        Amp = 1e-4
        mass_range = [0.3,2]
        tol = 1e-2
        run_type = "black_hole_mass_search"
        print("run_type = ", run_type)
        sim.A = Amp
        sim.record = run_params + "/record_{}.dat".format(l)
        sim.mass_search(l=l, mass_range = mass_range , tol = tol)
    else:
        for j in range(len(input_data)):
            launch_sim(input_data[j])

else:
    if sim.search == True:
        run_type = "black_hole_mass_search"
        sim.Amp = 0.
        tol = 1e-3
        cluster = True
        data_search = [[0.3,0.31,0.34],
              [0.4,0.43,0.46],
              [0.5,0.53,0.57],
              [0.6, 0.65,0.69],
              [0.7, 0.78,0.82],
              [0.8,0.87,0.90],
              [0.9,0.9,1.1],
              [1.,1.1, 1.15]]
        def launch_search(arr):
            l = arr[0]
            mass_range = [arr[1],arr[2]]
            sim.record = run_params + "/record_{}.dat".format(l)
            # dx = sim.cl/sim.nx
            # ratio = 0.5
            # sim.exc_i = int((ratio/dx)*((4*sim.cl*sim.initial_mass)/(sim.cl + 4*sim.initial_mass)))
            sim.mass_search(l=l, mass_range = mass_range , tol = tol)

        #--------------------------------------------------------------------------
        if __name__ == '__main__':
            print("theory = ",theory)
            if cluster:
                pool_nums = len(data_search)
            else:
                if len(data_search) >=6:
                    pool_nums = 6
                else :
                    pool_nums = len(data_search)


            print("pool_nums = ", pool_nums)

            t_start = time.time()

            print("Searching mass...")
            print("run_type = {}".format(run_type))
            print("tol = ", tol)

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
            print("Data saved at:{}".format(sim.out_dir))

            t_start = time.time()
            print("Starting multiprocessing pool..")
            pool = Pool(pool_nums)
            result = pool.map_async(launch_sim, input_data)
            pool.close()

            while True:
                if not result.ready():
                    print('We\'re not done yet, %s tasks to go!' % result._number_left)
                    time.sleep(5)
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

        # print("Ending caffeinate...")
        # subprocess.call("killall caffeinate" , shell = True)
        # print("Done. Exit.")
#===============================================================================
