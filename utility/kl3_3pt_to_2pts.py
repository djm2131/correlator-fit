#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

#################################################################
## A little utility function to take an all-mode averaged 3pt.
## function data file and split it up by src-sink separation.
##
## David Murphy (djm2131@columbia.edu)
## 12/13/2016
#################################################################

if len(sys.argv) != 5:
  print("Usage: ./kl3_3pt_to_2pts.py <data_stem> <traj_start> <traj_end> <traj_inc>")
  exit(0)

data_stem = sys.argv[1]
traj_start = int(sys.argv[2])
traj_end = int(sys.argv[3])
traj_incr = int(sys.argv[4])

for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = data_stem + "." + str(traj)
  print("Parsing trajectory " + fs + "...")
  data = np.genfromtxt(fs)
  seps = np.unique(data[:,0])
  for sep in seps:
    data_this_sep = data[data[:,0]==sep,:]
    fs_this_sep = data_stem + "." + str(int(sep)) + "." + str(traj)
    f = open(fs_this_sep,'w')
    for t in xrange(0,data_this_sep.shape[0]):
      line = "{0:3d} {1:1.10e} {2:1.10e} {3:1.10e} {4:1.10e} {5:1.10e} {6:1.10e} {7:1.10e} {8:1.10e} {9:1.10e} {10:1.10e}".format(  \
          int(data_this_sep[t,1]), data_this_sep[t,2], data_this_sep[t,3], \
          data_this_sep[t,4], data_this_sep[t,5], data_this_sep[t,6], data_this_sep[t,7],  \
          data_this_sep[t,8], data_this_sep[t,9], data_this_sep[t,10], data_this_sep[t,11])
      print(line, file=f)
    f.close()

print("done.")
