#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path

def time_avg(data, T):
  data_tavg = np.zeros((T,3))
  for t in xrange(0,T):
    data_tavg[t,0] = t
    data_tavg[t,1] = np.mean(data[t:T*T:T,2])
    data_tavg[t,2] = np.mean(data[t:T*T:T,3])
  return data_tavg

###################################################################
## A little utility function to process AMA data
##
## David Murphy (djm2131@columbia.edu)
## 04/17/2017
###################################################################

if len(sys.argv) != 8:
  print("Usage: ./process_ama.py <exact_fstem> <sloppy_fstem> <out_fstem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

exact_fstem  = sys.argv[1]
sloppy_fstem = sys.argv[2]
out_fstem    = sys.argv[3]
traj_start   = int(sys.argv[4])
traj_end     = int(sys.argv[5])
traj_incr    = int(sys.argv[6])
T            = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  print("Processing trajectory {0}...".format(traj), end="")
  tmp_exact = np.genfromtxt(exact_fstem + "." + str(traj))
  tmp_sloppy = np.genfromtxt(sloppy_fstem + "." + str(traj))
  tmp_exact[:,2:] -= tmp_sloppy[:,2:]
  tmp_exact = time_avg(tmp_exact,T)
  tmp_sloppy = time_avg(tmp_sloppy,T)
  tmp_exact[:,1:] += tmp_sloppy[:,1:]
  f = open(out_fstem + "." + str(traj), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.16e} {2:1.16e}".format(int(tmp_exact[t,0]),tmp_exact[t,1],tmp_exact[t,2])
    print(line, file=f)
  f.close()
  print("done.")

print("\nFinished!")
