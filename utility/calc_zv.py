#!/usr/bin/python3

import sys
import numpy as np
import os.path


###################################################################
## A little utility function to compute ZV from fit parameters.
##
## David Murphy (djm2131@columbia.edu)
## 06/28/2017
###################################################################

if len(sys.argv) != 4:
  print("Usage: ./calc_zv.py <ZVpsWW_stem> <ZpsWW_stem> <out_stem>")
  exit(0)

ZVpsWW_stem = sys.argv[1]
ZpsWW_stem  = sys.argv[2]
out_stem    = sys.argv[3]

# Fetch data
ZVpsWW_cv    = np.genfromtxt( ZVpsWW_stem + ".dat"       )
ZpsWW_cv     = np.genfromtxt( ZpsWW_stem  + ".dat"       )
ZVpsWW_jacks = np.genfromtxt( ZVpsWW_stem + "_jacks.dat" )
ZpsWW_jacks  = np.genfromtxt( ZpsWW_stem  + "_jacks.dat" )

with open(out_stem+".dat",'w') as f:
  ZV_cv = (ZpsWW_cv/ZVpsWW_cv)**2
  line = "{0:1.8e}".format(ZV_cv)
  print(line, file=f)

ZV_jacks = np.zeros((len(ZVpsWW_jacks)))
with open(out_stem+"_jacks.dat",'w') as f:
  for i in range(0,len(ZVpsWW_jacks)):
    ZV_jacks[i] = (ZpsWW_jacks[i]/ZVpsWW_jacks[i])**2
    line = "{0:1.8e}".format(ZV_jacks[i])
    print(line, file=f)

print("ZVps = {0} +/- {1}\n".format(ZV_cv,np.sqrt(len(ZV_jacks)-1.0)*np.std(ZV_jacks,ddof=0)))
