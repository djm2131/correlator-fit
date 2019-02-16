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
  print("Usage: ./process_2pion_ratio.py <two_pion_corr_stem> <out_stem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

Cpp_stem = sys.argv[1]
out_stem = sys.argv[2]
traj_start = int(sys.argv[3])
traj_end = int(sys.argv[4])
traj_incr = int(sys.argv[5])
T = int(sys.argv[6])

Ntraj = (traj_end - traj_start)/traj_incr + 1
Cpp = np.zeros((Ntraj,T))

# Fetch the two pion correlator data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = Cpp_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  Cpp[(traj-traj_start)/traj_incr,:] = data[:,1]

print("\nComputing jackknifed ratios...")

# Compute central values and error
R = np.zeros((T-2))
Cpp_jack = np.mean(Cpp, axis=0)
for t in xrange(1,T-1):
  R[t-1] = ( Cpp_jack[t+1] - Cpp_jack[t] ) / \
            ( Cpp_jack[t] - Cpp_jack[t-1] )
jacks = np.zeros((Ntraj,T-2))
for sample in xrange(0,Ntraj):
  Cpp_jack = np.mean( np.delete(Cpp, sample, axis=0), axis=0)
  for t in xrange(1,T-1):
    jacks[sample,t-1] = ( Cpp_jack[t+1] - Cpp_jack[t] ) / \
                        ( Cpp_jack[t] - Cpp_jack[t-1] )
weights = np.sqrt(Ntraj-1.0)*np.std(jacks, axis=0, ddof=0)

# Write central values to disk
f = open(out_stem, 'w')
for t in xrange(1,T-1):
  line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, R[t-1], weights[t-1])
  print(line, file=f)
f.close()

# Jackknife the data and write samples to disk
for sample in xrange(0,Ntraj):
  this_traj = traj_start + sample*traj_incr
  R = np.zeros((T-2))
  Cpp_jack = np.mean( np.delete(Cpp, sample, axis=0), axis=0)
  for t in xrange(1,T-1):
    R[t-1] = ( Cpp_jack[t+1] - Cpp_jack[t] ) / \
                ( Cpp_jack[t] - Cpp_jack[t-1] )
  jacks = np.zeros((Ntraj-1,T-2))
  for subsample in xrange(0,Ntraj-1):
    Cpp_jack_jack = np.mean( np.delete( np.delete(Cpp, sample, axis=0), subsample, axis=0), axis=0)
    for t in xrange(1,T-1):
      jacks[subsample,t-1] = ( Cpp_jack_jack[t+1] - Cpp_jack_jack[t] ) / \
                                ( Cpp_jack_jack[t] - Cpp_jack_jack[t-1] )
  weights = np.sqrt(Ntraj-2.0)*np.std(jacks, axis=0, ddof=0)
  f = open(out_stem + "." + str(this_traj), 'w')
  for t in xrange(1,T-1):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, R[t-1], weights[t-1])
    print(line, file=f)
  f.close()

print("done.")
