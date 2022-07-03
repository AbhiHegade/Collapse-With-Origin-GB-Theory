#!/usr/bin/env python
# coding: utf-8
import subprocess
id = "ah30"
subprocess.call("ssh {}@cc-login.campuscluster.illinois.edu".format(id),shell = True)
