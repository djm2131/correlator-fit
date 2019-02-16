#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

###################################################################
## A little utility function to process the raw Omega baryon
## correlator data
##
## David Murphy (djm2131@columbia.edu)
## 01/02/2017
###################################################################

## This function projects the omega two-point function onto its positive parity component.
## Used by import_omega.
def tr_proj_p(sss, trange):
  return 0.5*( sss[trange,1] + sss[trange,5] + sss[trange,11] + sss[trange,15] + \
               sss[trange,17] + sss[trange,21] + sss[trange,27] + sss[trange,31] )

## This function projects the omega two-point function onto its positive parity component.
## Used by import_omega.
def tr_proj_m(sss, trange):
  return 0.5*( sss[trange,1] - sss[trange,5] + sss[trange,11] - sss[trange,15] - \
               sss[trange,17] + sss[trange,21] - sss[trange,27] + sss[trange,31] )

## Computes negative parity correlator and then applies the transformation C(t) --> -C(Nt-t)
## Used by import_omega to average forward and backward propagating omega states
def get_backward_m_corr(dat, Nt):
  m_corr = tr_proj_m(dat,xrange(0,Nt))
  backward_m_corr = np.zeros((m_corr.shape[0]))
  backward_m_corr[0] = m_corr[0]
  for t in xrange(1,m_corr.shape[0]):
    backward_m_corr[t] = -m_corr[Nt-t]
  return backward_m_corr


if len(sys.argv) != 9:
  print("Usage: ./process_omega.py <stem_x_in> <stem_y_in> <stem_z_in> <stem_out> <traj_start> <traj_end> <traj_inc> <T>")
  exit(0)

stem_x_in  = sys.argv[1]
stem_y_in  = sys.argv[2]
stem_z_in  = sys.argv[3]
stem_out   = sys.argv[4]
traj_start = int(sys.argv[5])
traj_end   = int(sys.argv[6])
traj_incr  = int(sys.argv[7])
T          = int(sys.argv[8])

Ntraj = (traj_end - traj_start)/traj_incr + 1

# Fetch the data
for traj in xrange(traj_start, traj_end+1, traj_incr):

  print("Parsing trajectory " + str(traj) + "...")

  # Fetch data
  tmp_data = np.zeros((T,65))
  data_x = np.genfromtxt(stem_x_in + "." + str(traj))
  data_y = np.genfromtxt(stem_y_in + "." + str(traj))
  data_z = np.genfromtxt(stem_z_in + "." + str(traj))

  # Compute omega correlator
  for j in xrange(1,33):
    tmp_data[:,j] = 1.0/3.0*( data_x[:,j] + 2.0*data_x[:,j+32] + \
        data_y[:,j] + 2.0*data_y[:,j+32] + data_z[:,j] + 2.0*data_z[:,j+32] )
  back_corr = get_backward_m_corr(tmp_data, T)
  omega_corr = 0.5 * ( tr_proj_p(tmp_data, xrange(0,T)) + back_corr )

  # Write Omega correlator to disk
  fs = stem_out + "." + str(traj)
  f = open(fs,'w')
  for t in xrange(0,T):
    line = "{0:3d} {1:1.10e} {2:1.10e}".format(t, omega_corr[t], 0.0)
    print(line, file=f)
  f.close()

print("done.")
