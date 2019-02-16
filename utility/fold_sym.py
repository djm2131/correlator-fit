#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path

#################################################################
## A little utility function to fold a symmetric correlator
## about its midpoint.
##
## David Murphy (djm2131@columbia.edu)
## 12/28/2016
#################################################################

if len(sys.argv) != 7:
  print("Usage: ./fold_sym.py <in_stem> <out_stem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

in_stem = sys.argv[1]
out_stem = sys.argv[2]
traj_start = int(sys.argv[3])
traj_end = int(sys.argv[4])
traj_incr = int(sys.argv[5])
T = int(sys.argv[6])

if os.path.isfile(in_stem):
  print("Parsing trajectory " + in_stem + "...")
  data = np.genfromtxt(in_stem)
  if data.shape[0] != T:
    print("Data {0:s} does not match T={1:d}", in_stem, T)
    exit(-1)
  data_fold = np.zeros((T/2+1,3))
  data_fold[0,:] = data[0,:]
  for t in xrange(1,T/2+1):
    data_fold[t,0] = data[t,0]
    data_fold[t,1] = 0.5*(data[t,1]+data[T-t,1])
    data_fold[t,2] = 0.5*(data[t,2]+data[T-t,2])
  fout = out_stem
  f = open(fout,'w')
  for row in xrange(0,data_fold.shape[0]):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(int(data_fold[row,0]), data_fold[row,1], data_fold[row,2])
    print(line, file=f)
  f.close()

for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = in_stem + "." + str(traj)
  print("Parsing trajectory " + fs + "...")
  data = np.genfromtxt(fs)
  if data.shape[0] != T:
    print("Data {0:s} does not match T={1:d}", fs, T)
    exit(-1)
  data_fold = np.zeros((T/2+1,3))
  data_fold[0,:] = data[0,:]
  for t in xrange(1,T/2+1):
    data_fold[t,0] = data[t,0]
    data_fold[t,1] = 0.5*(data[t,1]+data[T-t,1])
    data_fold[t,2] = 0.5*(data[t,2]+data[T-t,2])
  fout = out_stem + "." + str(traj)
  f = open(fout,'w')
  for row in xrange(0,data_fold.shape[0]):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(int(data_fold[row,0]), data_fold[row,1], data_fold[row,2])
    print(line, file=f)
  f.close()

print("done.")
