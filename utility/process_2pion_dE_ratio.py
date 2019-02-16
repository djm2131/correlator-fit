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
  print("Usage: ./process_2pion_ratio.py <two_pion_corr_stem> <pion_corr_stem> <out_stem> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

two_pion_corr_stem = sys.argv[1]
pion_corr_stem = sys.argv[2]
out_stem = sys.argv[3]
traj_start = int(sys.argv[4])
traj_end = int(sys.argv[5])
traj_incr = int(sys.argv[6])
T = int(sys.argv[7])

Ntraj = (traj_end - traj_start)/traj_incr + 1
two_pion_corr = np.zeros((Ntraj,T))
pion_corr = np.zeros((Ntraj,T))

# Fetch the two pion correlator data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = two_pion_corr_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  two_pion_corr[(traj-traj_start)/traj_incr,:] = data[:,1]

# Fetch the pion correlator data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  fs = pion_corr_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  pion_corr[(traj-traj_start)/traj_incr,:] = data[:,1]

print("\nComputing jackknifed ratios...")

# Compute central values and error
R = np.zeros((T-1))
two_pion_jack = np.mean(two_pion_corr, axis=0)
pion_jack = np.mean(pion_corr, axis=0)
for t in xrange(0,T-1):
  R[t] = ( two_pion_jack[t] - two_pion_jack[t+1] ) / \
          ( pion_jack[t]**2 - pion_jack[t+1]**2 )
jacks = np.zeros((Ntraj,T-1))
for sample in xrange(0,Ntraj):
  two_pion_jack = np.mean( np.delete(two_pion_corr, sample, axis=0), axis=0)
  pion_jack = np.mean( np.delete(pion_corr, sample, axis=0), axis=0)
  for t in xrange(0,T-1):
    jacks[sample,t] = ( two_pion_jack[t] - two_pion_jack[t+1] ) / \
                        ( pion_jack[t]**2 - pion_jack[t+1]**2 )
weights = np.sqrt(Ntraj-1.0)*np.std(jacks, axis=0, ddof=0)

# Write central values to disk
f = open(out_stem, 'w')
for t in xrange(0,T-1):
  line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, R[t], weights[t])
  print(line, file=f)
f.close()

# Jackknife the data and write samples to disk
for sample in xrange(0,Ntraj):
  this_traj = traj_start + sample*traj_incr
  R = np.zeros((T-1))
  two_pion_jack = np.mean( np.delete(two_pion_corr, sample, axis=0), axis=0)
  pion_jack = np.mean( np.delete(pion_corr, sample, axis=0), axis=0)
  for t in xrange(0,T-1):
    R[t] = ( two_pion_jack[t] - two_pion_jack[t+1] ) / \
            ( pion_jack[t]**2 - pion_jack[t+1]**2 )
  jacks = np.zeros((Ntraj-1,T-1))
  for subsample in xrange(0,Ntraj-1):
    two_pion_jack_jack = np.mean( np.delete( np.delete(two_pion_corr, sample, axis=0), subsample, axis=0), axis=0)
    pion_jack_jack = np.mean( np.delete( np.delete(pion_corr, sample, axis=0), subsample, axis=0), axis=0)
    for t in xrange(0,T-1):
      jacks[subsample,t] = ( two_pion_jack_jack[t] - two_pion_jack_jack[t+1] ) / \
                            ( pion_jack_jack[t]**2 - pion_jack_jack[t+1]**2 )
  weights = np.sqrt(Ntraj-2.0)*np.std(jacks, axis=0, ddof=0)
  f = open(out_stem + "." + str(this_traj), 'w')
  for t in xrange(0,T-1):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, R[t], weights[t])
    print(line, file=f)
  f.close()

print("done.")
