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

if len(sys.argv) != 8:
  print("Usage: ./process_zv_corr.py <zv_stem_in> <zv_ref_stem_in> <zv_stem_out> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

zv_stem_in     = sys.argv[1]
zv_ref_stem_in = sys.argv[2]
zv_stem_out    = sys.argv[3]
traj_start     = int(sys.argv[4])
traj_end       = int(sys.argv[5])
traj_incr      = int(sys.argv[6])
T              = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
ZV    = np.zeros((Ntraj,T,T))

for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  print("Parsing trajectory " + str(traj) + ":")
  for sep in xrange(0,T):
    fs = zv_stem_in + "." + str(sep) + "." + str(traj)
    fs_ref = zv_ref_stem_in + "." + str(sep) + "." + str(traj)
    if os.path.isfile(fs):
      print("\tSource-sink separation " + str(sep) + "...")
      data = np.genfromtxt(fs)
      data_ref = np.genfromtxt(fs_ref)
      ZV[idx,:,sep] = 0.5*( -data[:,7] + data_ref[:,7] )
      f = open(zv_stem_out + "." + str(sep) + "." + str(traj), 'w')
      for t in xrange(0,T):
        line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, ZV[idx,t,sep], 0.0)
        print(line, file=f)
      f.close()
  print("")

print("done.")
