#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np
import os.path


###################################################################
## A little utility function to process raw Z_{V} correlator data
## and compute each jackknife sample of ratio we conventially fit.
##
## David Murphy (djm2131@columbia.edu)
## 10/05/2016
###################################################################

if len(sys.argv) != 10:
  print("Usage: ./calc_zv_eff_corr.py <Z_stem> <m_stem> <ZV_stem> <out_stem> <sep_start> <sep_end> <sep_inc> <L> <T>")
  exit(0)

Z_stem     = sys.argv[1]
m_stem     = sys.argv[2]
ZV_stem    = sys.argv[3]
out_stem   = sys.argv[4]
sep_start  = int(sys.argv[5])
sep_end    = int(sys.argv[6])
sep_inc    = int(sys.argv[7])
L          = int(sys.argv[8])
T          = int(sys.argv[9])

# Fetch data
Z_cv     = np.genfromtxt( Z_stem  + ".dat"       )
m_cv     = np.genfromtxt( m_stem  + ".dat"       )
ZV_cv    = np.genfromtxt( ZV_stem + ".dat"       )
Z_jacks  = np.genfromtxt( Z_stem  + "_jacks.dat" )
m_jacks  = np.genfromtxt( m_stem  + "_jacks.dat" )
ZV_jacks = np.genfromtxt( ZV_stem + "_jacks.dat" )

V     = L**3
Ntraj = Z_jacks.shape[0]

ZV_corr_cv    = 0.0
ZV_corr_jacks = np.zeros((Ntraj))
for sep in xrange(sep_start, sep_end+1, sep_inc):

  ZV_corr_cv = 0.5*Z_cv**2/m_cv/V/ZV_cv * np.exp(-m_cv*sep)
  for jack in xrange(0,Ntraj):
    ZV_corr_jacks[jack] = 0.5*Z_jacks[jack]**2/m_jacks[jack]/V/ZV_jacks[jack] * np.exp(-m_jacks[jack]*sep)

  # Write central value to disk
  f_cv = open(out_stem + "." + str(sep) + ".dat", 'w')
  print("{0:1.10e}".format(ZV_corr_cv), file=f_cv)
  f_cv.close()

  # Write jackknife samples to disk
  f_jacks = open(out_stem + "_jacks." + str(sep) + ".dat", 'w')
  for jack in xrange(0,Ntraj):
    print("{0:1.10e}".format(ZV_corr_jacks[jack]), file=f_jacks)
  f_jacks.close()
