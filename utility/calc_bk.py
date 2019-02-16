#!/usr/bin/python3

import sys
import numpy as np
import os.path


###################################################################
## A little utility function to compute pseudoscalar decay
## constants from fit parameters.
##
## David Murphy (djm2131@columbia.edu)
## 06/28/2017
###################################################################

if len(sys.argv) != 7:
  print("Usage: ./calc_bk.py <mK_stem> <ZKWW_stem> <ZBKWW_stem> <ZA_stem> <fK_stem> <out_stem>")
  exit(0)

mK_stem    = sys.argv[1]
ZKWW_stem  = sys.argv[2]
ZBKWW_stem = sys.argv[3]
ZA_stem    = sys.argv[4]
fK_stem    = sys.argv[5]
out_stem   = sys.argv[6]

# Fetch data
mK_cv       = np.genfromtxt( mK_stem    + ".dat"       )
ZKWW_cv     = np.genfromtxt( ZKWW_stem  + ".dat"       )
ZBKWW_cv    = np.genfromtxt( ZBKWW_stem + ".dat"       )
ZA_cv       = np.genfromtxt( ZA_stem    + ".dat"       )
fK_cv       = np.genfromtxt( fK_stem    + ".dat"       )
mK_jacks    = np.genfromtxt( mK_stem    + "_jacks.dat" )
ZKWW_jacks  = np.genfromtxt( ZKWW_stem  + "_jacks.dat" )
ZBKWW_jacks = np.genfromtxt( ZBKWW_stem + "_jacks.dat" )
ZA_jacks    = np.genfromtxt( ZA_stem    + "_jacks.dat" )
fK_jacks    = np.genfromtxt( fK_stem    + "_jacks.dat" )

with open(out_stem+".dat",'w') as f:
  BK_cv = 3.0/(4.0*mK_cv) * (ZBKWW_cv*ZA_cv/ZKWW_cv/fK_cv)**2
  line = "{0:1.8e}".format(BK_cv)
  print(line, file=f)

BK_jacks = np.zeros((len(mK_jacks)))
with open(out_stem+"_jacks.dat",'w') as f:
  for i in range(0,len(mK_jacks)):
    BK_jacks[i] = 3.0/(4.0*mK_jacks[i]) * (ZBKWW_jacks[i]*ZA_jacks[i]/ZKWW_jacks[i]/fK_jacks[i])**2
    line = "{0:1.8e}".format(BK_jacks[i])
    print(line, file=f)

print("B_K = {0} +/- {1}\n".format(BK_cv,np.sqrt(len(BK_jacks)-1.0)*np.std(BK_jacks,ddof=0)))
