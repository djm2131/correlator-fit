#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path


###################################################################
## A little utility function to process raw B_{K} correlator data
## and compute each jackknife sample of ratio we conventially fit.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
###################################################################

if len(sys.argv) != 9:
  print("Usage: ./calc_bk_eff_corr.py <ZBKWW_stem> <mK_stem> <out_stem> <sep_start> <sep_end> <sep_inc> <L> <T>")
  exit(0)

ZBKWW_stem = sys.argv[1]
mK_stem    = sys.argv[2]
out_stem   = sys.argv[3]
sep_start  = int(sys.argv[4])
sep_end    = int(sys.argv[5])
sep_inc    = int(sys.argv[6])
L          = int(sys.argv[7])
T          = int(sys.argv[8])

# Fetch data
ZBKWW_cv    = np.genfromtxt( ZBKWW_stem + ".dat"       )
mK_cv       = np.genfromtxt( mK_stem    + ".dat"       )
ZBKWW_jacks = np.genfromtxt( ZBKWW_stem + "_jacks.dat" )
mK_jacks    = np.genfromtxt( mK_stem    + "_jacks.dat" )
V           = L**3
Ntraj       = mK_jacks.shape[0]

BK_corr_cv    = 0.0
BK_corr_jacks = np.zeros((Ntraj))
for sep in xrange(sep_start, sep_end+1, sep_inc):

  BK_corr_cv = ZBKWW_cv**2 / (2.0*mK_cv*V) * np.exp(-mK_cv*sep)
  for jack in xrange(0,Ntraj):
    BK_corr_jacks[jack] = ZBKWW_jacks[jack]**2 / (2.0*mK_jacks[jack]*V) * np.exp(-mK_jacks[jack]*sep)

  # Write central value to disk
  f_cv = open(out_stem + "." + str(sep) + ".dat", 'w')
  print("{0:1.10e}".format(BK_corr_cv), file=f_cv)
  f_cv.close()

  # Write jackknife samples to disk
  f_jacks = open(out_stem + "_jacks." + str(sep) + ".dat", 'w')
  for jack in xrange(0,Ntraj):
    print("{0:1.10e}".format(BK_corr_jacks[jack]), file=f_jacks)
  f_jacks.close()
