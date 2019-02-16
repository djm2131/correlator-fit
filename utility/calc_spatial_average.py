#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

#################################################################
## A little utility function to fold a symmetric correlator
## about its midpoint.
##
## David Murphy (djm2131@columbia.edu)
## 12/28/2016
#################################################################

if len(sys.argv) != 8:
  print("Usage: ./spatial_avg.py <in_stem_x> <in_stem_y> <in_stem_z> <out_stem> <traj_start> <traj_end> <traj_inc>")
  exit(0)

in_stem_x = sys.argv[1]
in_stem_y = sys.argv[2]
in_stem_z = sys.argv[3]
out_stem = sys.argv[4]
traj_start = int(sys.argv[5])
traj_end = int(sys.argv[6])
traj_incr = int(sys.argv[7])

for traj in xrange(traj_start, traj_end+1, traj_incr):
  print("Parsing trajectory " + str(traj) + "...")
  fs_x = in_stem_x + "." + str(traj)
  fs_y = in_stem_y + "." + str(traj)
  fs_z = in_stem_z + "." + str(traj)
  data_x = np.genfromtxt(fs_x)
  data_y = np.genfromtxt(fs_y)
  data_z = np.genfromtxt(fs_z)
  T = data_x.shape[0]
  data_avg = np.zeros((T,3))
  data_avg[:,0] = data_x[:,0]
  data_avg[:,1] = ( data_x[:,1] + data_y[:,1] + data_z[:,1] ) / 3.0
  data_avg[:,2] = ( data_x[:,2] + data_y[:,2] + data_z[:,2] ) / 3.0
  fout = out_stem + ".avg." + str(traj)
  f = open(fout,'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(int(data_avg[t,0]), data_avg[t,1], data_avg[t,2])
    print(line, file=f)
  f.close()

print("done.")
