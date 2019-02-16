#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path


###################################################################
## A little utility function to process raw Z_{V} correlator data
## and compute each jackknife sample of ratio we conventially fit.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
###################################################################

if len(sys.argv) != 7:
  print("Usage: ./process_bk.py <bk_stem_in> <bk_stem_out> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

bk_stem_in  = sys.argv[1]
bk_stem_out = sys.argv[2]
traj_start  = int(sys.argv[3])
traj_end    = int(sys.argv[4])
traj_incr   = int(sys.argv[5])
T           = int(sys.argv[6])

Ntraj = (traj_end - traj_start)/traj_incr + 1
bk    = np.zeros((Ntraj,T,T))

for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  print("Parsing trajectory " + str(traj) + ":")
  for sep in xrange(0,T):
    fs = bk_stem_in + "." + str(sep) + "." + str(traj)
    if os.path.isfile(fs):
      print("\tSource-sink separation " + str(sep) + "...")
      data = np.genfromtxt(fs)
      bk[idx,:,sep] = -2.0*( data[:,1] - data[:,3] )
      f = open(bk_stem_out + "." + str(sep) + "." + str(traj), 'w')
      for t in xrange(0,T):
        line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, bk[idx,t,sep], 0.0)
        print(line, file=f)
      f.close()
  print("")

print("done.")
