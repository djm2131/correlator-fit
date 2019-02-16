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
  print("Usage: ./calc_ps_decay.py <mps_stem> <ZpsWW_stem> <ZpsAP_stem> <ZA_stem> <L> <out_stem>")
  exit(0)

mps_stem   = sys.argv[1]
ZpsWW_stem = sys.argv[2]
ZpsAP_stem = sys.argv[3]
ZA_stem   = sys.argv[4]
L         = float(sys.argv[5])
out_stem  = sys.argv[6]

# Fetch data
mps_cv      = np.genfromtxt( mps_stem   + ".dat"       )
ZpsWW_cv    = np.genfromtxt( ZpsWW_stem + ".dat"       )
ZpsAP_cv    = np.genfromtxt( ZpsAP_stem + ".dat"       )
ZA_cv      = np.genfromtxt( ZA_stem   + ".dat"       )
mps_jacks   = np.genfromtxt( mps_stem   + "_jacks.dat" )
ZpsWW_jacks = np.genfromtxt( ZpsWW_stem + "_jacks.dat" )
ZpsAP_jacks = np.genfromtxt( ZpsAP_stem + "_jacks.dat" )
ZA_jacks   = np.genfromtxt( ZA_stem   + "_jacks.dat" )

with open(out_stem+".dat",'w') as f:
  fps_cv = (ZA_cv*ZpsAP_cv**2) / (mps_cv*L**3*ZpsWW_cv)
  line = "{0:1.8e}".format(fps_cv)
  print(line, file=f)

fps_jacks = np.zeros((len(mps_jacks)))
with open(out_stem+"_jacks.dat",'w') as f:
  for i in range(0,len(mps_jacks)):
    fps_jacks[i] = (ZA_jacks[i]*ZpsAP_jacks[i]**2) / (mps_jacks[i]*L**3*ZpsWW_jacks[i])
    line = "{0:1.8e}".format(fps_jacks[i])
    print(line, file=f)

print("fps = {0} +/- {1}\n".format(fps_cv,np.sqrt(len(fps_jacks)-1.0)*np.std(fps_jacks,ddof=0)))
