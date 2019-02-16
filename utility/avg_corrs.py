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
  print("Usage: ./process_tbc.py <stem_1> <stem_2> <stem_out> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

stem_1     = sys.argv[1]
stem_2     = sys.argv[2]
stem_out   = sys.argv[3]
traj_start = int(sys.argv[4])
traj_end   = int(sys.argv[5])
traj_incr  = int(sys.argv[6])
T          = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  print("Processing trajectory {0}...".format(traj), end="")
  C1 = np.genfromtxt(stem_1 + "." + str(traj))
  C2 = np.genfromtxt(stem_2 + "." + str(traj))
  f = open(stem_out + "." + str(traj), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.16e} {2:1.16e}".format(t,0.5*(C1[t,1]+C2[t,1]),0.5*(C1[t,2]+C2[t,2]))
    print(line, file=f)
  f.close()
  print("done.")

print("\nFinished!")
