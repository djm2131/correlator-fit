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
num = np.zeros((Ntraj,T))
denom = np.zeros((Ntraj,T))

# Fetch the data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = in_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  num[(traj-traj_start)/traj_incr,:] = data[:,3]
  denom[(traj-traj_start)/traj_incr,:] = data[:,1]

print("\nComputing jackknifed ratios...")

# Compute central values and error
mres = np.mean(num, axis=0) / np.mean(denom, axis=0)
jacks = np.zeros((Ntraj,T))
num_jack = np.zeros((Ntraj-1,T))
denom_jack = np.zeros((Ntraj-1,T))
for sample in xrange(0,Ntraj):
  num_jack = np.delete(num, sample, axis=0)
  denom_jack = np.delete(denom, sample, axis=0)
  jacks[sample,:] = np.mean(num_jack, axis=0) / np.mean(denom_jack, axis=0)
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
  num_jack = np.delete(num, sample, axis=0)
  denom_jack = np.delete(denom, sample, axis=0)
  mres = np.mean(num_jack, axis=0) / np.mean(denom_jack, axis=0)
  jacks = np.zeros((Ntraj-1,T))
  for subsample in xrange(0,Ntraj-1):
    num_jack_jack = np.delete(num_jack, subsample, axis=0)
    denom_jack_jack = np.delete(denom_jack, subsample, axis=0)
    jacks[subsample,:] = np.mean(num_jack_jack, axis=0) / np.mean(denom_jack_jack, axis=0)
  weights = np.sqrt(Ntraj-2.0)*np.std(jacks, axis=0, ddof=0)
  f = open(out_stem + "." + str(this_traj), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, mres[t], weights[t])
    print(line, file=f)
  f.close()

print("done.")
