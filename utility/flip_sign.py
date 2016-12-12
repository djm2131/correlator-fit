#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

#################################################################
## A little utility function to multiply correlator data by -1.
## Useful for fixing sign errors in the data generation.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
#################################################################

if len(sys.argv) != 5:
  print("Usage: ./flip_sign.py <data_stem> <traj_start> <traj_end> <traj_inc>")
  exit(0)

data_stem = sys.argv[1]
traj_start = int(sys.argv[2])
traj_end = int(sys.argv[3])
traj_incr = int(sys.argv[4])

for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = data_stem + "." + str(traj)
  print("Parsing trajectory " + fs + "...")
  data = np.genfromtxt(fs)
  f = open(fs,'w')
  for row in xrange(0,data.shape[0]):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(int(data[row,0]), -1.0*data[row,1], -1.0*data[row,2])
    print(line, file=f)
  f.close()

print("done.")
