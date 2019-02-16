#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path

###################################################################
## A little utility function to process AMA data
##
## David Murphy (djm2131@columbia.edu)
## 04/17/2017
###################################################################

if len(sys.argv) != 8:
  print("Usage: ./process_tbc.py <FF_fstem> <BB_fstem> <out_fstem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

FF_fstem   = sys.argv[1]
BB_fstem   = sys.argv[2]
out_fstem  = sys.argv[3]
traj_start = int(sys.argv[4])
traj_end   = int(sys.argv[5])
traj_incr  = int(sys.argv[6])
T          = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  print("Processing trajectory {0}...".format(traj), end="")
  FF = np.genfromtxt(FF_fstem + "." + str(traj))
  BB = np.genfromtxt(BB_fstem + "." + str(traj))
  f = open(out_fstem + "." + str(traj), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.16e} {2:1.16e}".format(t,FF[t,1],FF[t,2])
    print(line, file=f)
  for t in xrange(0,T):
    line = "{0:3d} {1:1.16e} {2:1.16e}".format(t+T,BB[t,1],BB[t,2])
    print(line, file=f)
  f.close()
  print("done.")

print("\nFinished!")
