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

if len(sys.argv) != 12:
  print("Usage: ./process_zv_ratio.py <pion_stem> <zv_stem_in> <zv_ref_stem_in> <zv_stem_out> <traj_start> <traj_end> <traj_inc> <L> <T> <Z_file> <m_file>")
  exit(0)

pion_stem      = sys.argv[1]
zv_stem_in     = sys.argv[2]
zv_stem_ref_in = sys.argv[3]
zv_stem_out    = sys.argv[4]
traj_start     = int(sys.argv[5])
traj_end       = int(sys.argv[6])
traj_incr      = int(sys.argv[7])
L              = int(sys.argv[8])
T              = int(sys.argv[9])
Z_file         = sys.argv[10]
m_file         = sys.argv[11]

Ntraj    = (traj_end - traj_start)/traj_incr + 1
pion     = np.zeros((Ntraj,T))
ZV       = np.zeros((Ntraj,T,T))

# Fetch fit parameter jackknife samples
V = L**3
Z = np.genfromtxt(Z_file)
m = np.genfromtxt(m_file)

# Fetch the pion data and subtract thermal state
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  fs = pion_stem + "." + str(traj)
  print("Parsing file " + fs + "...")
  data = np.genfromtxt(fs)
  pion[idx,:] = data[:,1]
  for t in xrange(0,T):
    # print(Z[idx]**2/(2.0*m[idx]*V)*np.exp(-m[idx]*(T-t)))
    # print(0.5*pion[idx,T/2]*np.exp(m[idx]*(t-T/2)))
    # pion[idx,t] -= Z[idx]**2/(2.0*m[idx]*V)*np.exp(-m[idx]*(T-t))
    pion[idx,t] -= 0.5*pion[idx,T/2]*np.exp(m[idx]*(t-T/2))
print("done.\n")

# Fetch the ZV data
for traj in xrange(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  for sep in xrange(0,T):
    fs = zv_stem_in + "." + str(sep) + "." + str(traj)
    fs_ref = zv_stem_ref_in + "." + str(sep) + "." + str(traj)
    if os.path.isfile(fs):
      print("Parsing file " + fs + "...")
      data = np.genfromtxt(fs)
      data_ref = np.genfromtxt(fs_ref)
      ZV[idx,:,sep] = 0.5*( -data[:,7] + data_ref[:,7] )

print("\nComputing jackknifed ratios...")

# Compute central value and error
ZV_ratio = np.zeros((T))
num = np.mean(pion, axis=0)
for sep in xrange(0,T):

  # Skip separations with no data
  if ZV[0,1,sep] == 0.0:
    continue

  # Compute ZV ratio
  den = np.mean(ZV[:,:,sep], axis=0)
  for t in xrange(0,T):
    ZV_ratio[t] = num[sep] / den[t]
  jacks = np.zeros((Ntraj,T))
  for sample in xrange(0,Ntraj):
    num = np.mean( np.delete(pion, sample, axis=0), axis=0)
    den = np.mean( np.delete(ZV[:,:,sep], sample, axis=0), axis=0)
    for t in xrange(0,T):
      jacks[sample,t] = num[sep] / den[t]
  weights = np.sqrt(Ntraj-1.0)*np.std(jacks, axis=0, ddof=0)

  # Write to disk
  f = open(zv_stem_out + "." + str(sep), 'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, ZV_ratio[t], weights[t])
    print(line, file=f)
  f.close()

# Jackknife the data and write samples to disk
for sep in xrange(0,T):

  # Skip separations with no data
  if ZV[0,1,sep] == 0.0:
    continue

  for sample in xrange(0,Ntraj):
    this_traj = traj_start + sample*traj_incr
    num_jacks = np.delete(pion, sample, axis=0)
    den_jacks = np.delete(ZV[:,:,sep], sample, axis=0)
    num = np.mean(num_jacks, axis=0)
    den = np.mean(den_jacks, axis=0)
    for t in xrange(0,T):
      ZV_ratio[t] = num[sep] / den[t]
    jacks = np.zeros((Ntraj-1,T))
    for subsample in xrange(0,Ntraj-1):
      num = np.mean( np.delete(num_jacks, subsample, axis=0), axis=0)
      den = np.mean( np.delete(den_jacks, subsample, axis=0), axis=0)
      for t in xrange(0,T):
        jacks[subsample,t] = num[sep] / den[t]
    weights = np.sqrt(Ntraj-2.0)*np.std(jacks, axis=0, ddof=0)
    f = open(zv_stem_out + "." + str(sep) + "." + str(this_traj), 'w')
    for t in xrange(0,T):
      line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, ZV_ratio[t], weights[t])
      print(line, file=f)
    f.close()

print("done.")
