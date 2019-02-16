#!/usr/bin/python3

import sys
import numpy as np
import os.path


###################################################################
## A little utility function to compute Omega baryon effective
## mass function for plotting.
##
## David Murphy (djm2131@columbia.edu)
## 01/14/2017
###################################################################

if len(sys.argv) != 4:
  print("Usage: ./calc_Epipi.py <mp_stem> <dEpp_stem> <out_stem>")
  exit(0)

mp_stem   = sys.argv[1]
dEpp_stem = sys.argv[2]
out_stem  = sys.argv[3]

# Fetch data
mp_cv      = np.genfromtxt( mp_stem   + ".dat"       )
dEpp_cv    = np.genfromtxt( dEpp_stem + ".dat"       )
mp_jacks   = np.genfromtxt( mp_stem   + "_jacks.dat" )
dEpp_jacks = np.genfromtxt( dEpp_stem + "_jacks.dat" )

with open(out_stem+".dat",'w') as f:
  line = "{0:1.8e}".format(dEpp_cv + 2.0*mp_cv)
  print(line, file=f)

with open(out_stem+"_jacks.dat",'w') as f:
  for i in range(0,len(dEpp_jacks)):
    line = "{0:1.8e}".format(dEpp_jacks[i] + 2.0*mp_jacks[i])
    print(line, file=f)
