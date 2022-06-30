#!/usr/bin/env python
# coding: utf-8
from datetime import datetime
import os, sys, time, shutil, subprocess

class Sim:
    def __init__(self):
        self.home_dir = str(os.getcwd())
#===============================================================
    def make_output_dir(self):
        current_time = datetime.now()
        self.output_dir = self.out_dir + "/" + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"+ str(current_time.day) +"_"+ str(current_time.hour) + "_" + str(current_time.hour)+"_"+str(current_time.minute)+"_"+str(current_time.second)
        self.output_dir = self.output_dir + "_A_" + "{:.2e}".format(self.A) +"_l_" + "{:.2e}".format(self.l)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
#===============================================================
    def write_sim_params(self):
        with open(self.output_dir+'/sim_params.txt','w') as f:
            attrs= vars(self)
            for param in attrs:
                f.write('{} {}\n'.format(param,attrs[param]))
#===============================================================
    def make_output_file(self):
        self.output_file= self.output_dir+'/'+'output.out'
        with open(self.output_file, 'w') as f:
            pass
#===============================================================
    def launch(self):
        self.make_output_dir()
        self.make_output_file()
        self.write_sim_params()
        subprocess.call('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
        # subprocess.Popen('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
