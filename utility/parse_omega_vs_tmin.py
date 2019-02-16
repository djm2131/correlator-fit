#!/usr/bin/python3

import numpy as np

box_min = 3
box_max = 6
wall_min = 5
wall_max = 8
mO = np.zeros((((box_max-box_min+1)*(wall_max-wall_min+1)),4))
mOp = np.zeros((((box_max-box_min+1)*(wall_max-wall_min+1)),4))

idx = 0
results_flag = False
for tmin_box in range(box_min, box_max+1):
  for tmin_wall in range(wall_min, wall_max+1):
    mO[idx,0] = tmin_box
    mOp[idx,0] = tmin_box
    mO[idx,1] = tmin_wall
    mOp[idx,1] = tmin_wall
    with open("logs/fit_omega_{0}_{1}.log".format(tmin_box, tmin_wall), 'r') as f:
      for line in f:
        if "Fit results" in line:
          results_flag = True
          continue
        if results_flag and "mO" in line and "mOp" not in line:
          mO[idx,2] = float( line.split()[-3] )
          mO[idx,3] = float( line.split()[-1] )
        if results_flag and "mOp" in line:
          mOp[idx,2] = float( line.split()[-3] )
          mOp[idx,3] = float( line.split()[-1] )
          results_flag = False
          idx += 1

np.savetxt("mO_vs_tmin.dat", mO)
np.savetxt("mOp_vs_tmin.dat", mOp)
