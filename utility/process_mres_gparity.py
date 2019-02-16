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

if len(sys.argv) != 8:
  print("Usage: ./process_mres.py <J5q_stem> <J5_stem> <out_stem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

J5q_stem = sys.argv[1]
J5_stem = sys.argv[2]
out_stem = sys.argv[3]
traj_start = int(sys.argv[4])
traj_end = int(sys.argv[5])
traj_incr = int(sys.argv[6])
T = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
J5q = np.zeros((Ntraj,T))
J5 = np.zeros((Ntraj,T))

# Fetch the data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = J5q_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  J5q[(traj-traj_start)/traj_incr,:] = data[:,1]
  fs = J5_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  J5[(traj-traj_start)/traj_incr,:] = data[:,1]

print("\nComputing jackknifed ratios...")

# Compute central values and error
mres = np.mean(J5q, axis=0) / np.mean(J5, axis=0)
jacks = np.zeros((Ntraj,T))
J5q_jack = np.zeros((Ntraj-1,T))
J5_jack = np.zeros((Ntraj-1,T))
for sample in xrange(0,Ntraj):
  J5q_jack = np.delete(J5q, sample, axis=0)
  J5_jack = np.delete(J5, sample, axis=0)
  jacks[sample,:] = np.mean(J5q_jack, axis=0) / np.mean(J5_jack, axis=0)
weights = np.sqrt(Ntraj-1.0)*np.std(jacks, axis=0, ddof=0)

# Write central values to disk
f = open(out_stem, 'w')
for t in xrange(0,T):
  line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, mres[t], weights[t])
  print(line, file=f)
f.close()

# Jackknife the data and write samples to disk
for sample in xrange(0,Ntraj):
  this_traj = traj_start + sample*traj_incr
  J5q_jack = np.delete(J5q, sample, axis=0)
  J5_jack = np.delete(J5, sample, axis=0)
  mres = np.mean(J5q_jack, axis=0) / np.mean(J5_jack, axis=0)
  jacks = np.zeros((Ntraj-1,T))
  for subsample in xrange(0,Ntraj-1):
    J5q_jack_jack = np.delete(J5q_jack, subsample, axis=0)
    J5_jack_jack = np.delete(J5_jack, subsample, axis=0)
    jacks[subsample,:] = np.mean(J5q_jack_jack, axis=0) / np.mean(J5_jack_jack, axis=0)
  weights = np.sqrt(Ntraj-2.0)*np.std(jacks, axis=0, ddof=0)
  f = open(out_stem + "." + str(this_traj), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, mres[t], weights[t])
    print(line, file=f)
  f.close()

print("done.")
