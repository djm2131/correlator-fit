#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path

def compute_za_ratio(za_num, za_denom):

  tt = len(za_num) - 2
  za_ratio = np.zeros((tt))

  for t in xrange(1,tt+1):
    za_ratio[t-1] = 0.5*( 0.5*(za_num[t]+za_num[t-1])/za_denom[t] + \
        2.0*za_num[t]/(za_denom[t+1]+za_denom[t]) )

  return za_ratio

###################################################################
## A little utility function to process raw Z_{V} correlator data
## and compute each jackknife sample of ratio we conventially fit.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
###################################################################

if len(sys.argv) != 9:
  print("Usage: ./process_bk_ratio.py <KAP_stem> <BK_stem_in> <BK_stem_out> <traj_start> <traj_end> <traj_inc> <L> <T>")
  exit(0)

KAP_stem    = sys.argv[1]
BK_stem_in  = sys.argv[2]
BK_stem_out = sys.argv[3]
traj_start  = int(sys.argv[4])
traj_end    = int(sys.argv[5])
traj_incr   = int(sys.argv[6])
L           = int(sys.argv[7])
T           = int(sys.argv[8])

Ntraj = (traj_end - traj_start)/traj_incr + 1
CKAP  = np.zeros((Ntraj,T))
BK    = np.zeros((Ntraj,T,T))

# Fetch the kaon <AP> data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  fs = KAP_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  CKAP[idx,:] = data[:,1]
print("done.\n")

# Fetch the BK data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  for sep in xrange(0,T):
    fs = BK_stem_in + "." + str(sep) + "." + str(traj)
    if os.path.isfile(fs):
      print("Parsing file " + fs + "...")
      data = np.genfromtxt(fs)
      BK[idx,:,sep] = 2.0*( data[:,1] - data[:,3] )

print("\nComputing jackknifed ratios...")

# Compute central value and error
BK_ratio = np.zeros((T))
std_err  = np.zeros((T))
for sep in xrange(0,T):

  # Skip separations with no data
  if (sep == 0) or (sep%2 != 0) or (BK[0,1,sep] == 0.0):
    BK_ratio[sep] = std_err[sep] = np.nan
    continue

  # Compute BK ratio
  num = np.mean(BK[:,:,sep], axis=0)
  den = np.mean(CKAP, axis=0)
  num = num[sep/2]
  for t in xrange(0,T):
    BK_ratio[sep] = 3.0/8.0 * L**3 * num / den[sep/2] / den[T-sep/2]
  jacks = np.zeros((Ntraj))
  for sample in xrange(0,Ntraj):
    num = np.mean( np.delete(BK[:,:,sep], sample, axis=0), axis=0)
    num = num[sep/2]
    den = np.mean( np.delete(CKAP, sample, axis=0), axis=0)
    jacks[sample] = 3.0/8.0 * L**3 * num / den[sep/2] / den[T-sep/2]
  std_err[sep] = np.sqrt(Ntraj-1.0)*np.std(jacks, axis=0, ddof=0)

# Write to disk
f = open(BK_stem_out, 'w')
for sep in xrange(0,T):
  if sep == 0:
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(sep, np.nan, np.nan)
  else:
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(sep, BK_ratio[sep], std_err[sep])
  print(line, file=f)
f.close()

print("done.")
