#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

###################################################################
## A little utility function to process raw mres correlator data
## and compute each jackknife sample of ratio we conventially fit.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
###################################################################

if len(sys.argv) != 7:
  print("Usage: ./process_mres.py <in_stem> <out_stem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

in_stem = sys.argv[1]
out_stem = sys.argv[2]
traj_start = int(sys.argv[3])
traj_end = int(sys.argv[4])
traj_incr = int(sys.argv[5])
T = int(sys.argv[6])

Ntraj = (traj_end - traj_start)/traj_incr + 1
D_diagram = np.zeros((Ntraj,T))
C_diagram = np.zeros((Ntraj,T))

# Fetch the data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = in_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  D_diagram[(traj-traj_start)/traj_incr,:] = data[:,1]
  C_diagram[(traj-traj_start)/traj_incr,:] = data[:,3]

# Write two pion correlator to disk
tt = 0
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = out_stem + "." + str(traj)
  f = open(fs,'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, 2.0*(D_diagram[tt,t]-C_diagram[tt,t]), 0.0)
    print(line, file=f)
  f.close()
  tt += 1

print("done.")
