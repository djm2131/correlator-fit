#!/usr/bin/python

from __future__ import print_function
import sys
import numpy as np

#################################################################
## A little utility function to bin correlation functions.
##
## David Murphy (djm2131@columbia.edu)
## 12/30/2016
#################################################################

if len(sys.argv) != 6:
  print("Usage: ./bin_3pt.py <data_stem> <traj_start> <traj_end> <traj_inc> <bin_size>")
  exit(0)

data_stem  = sys.argv[1]
traj_start = int(sys.argv[2])
traj_end   = int(sys.argv[3])
traj_incr  = int(sys.argv[4])
bin_size   = int(sys.argv[5])

# Get number of bins
Ntraj = (traj_end-traj_start)/traj_incr + 1
if Ntraj % bin_size != 0:
  print("Error: bin size {0:d} does not evenly divide trajectories.".format(bin_size))
Nbins = Ntraj / bin_size

# Get data file dimensions
fs = data_stem + "." + str(traj_start)
tmp = np.genfromtxt(fs)
Nrows = tmp.shape[0]
Ncols = tmp.shape[1]

data = np.zeros((Nrows,Ncols,Ntraj))
for traj in range(traj_start, traj_end+1, traj_incr):
  idx = (traj-traj_start)/traj_incr
  fs = data_stem + "." + str(traj)
  print("Parsing trajectory " + fs + "...")
  data[:,:,idx] = np.genfromtxt(fs)
print("done.\n")

data_binned = np.zeros((Nrows,Ncols,Nbins))
for traj in range(traj_start, traj_end+1, bin_size*traj_incr):
  print("Binning trajectories " + str(traj) + "-" + \
          str((bin_size-1)*traj_incr + traj) + "...")
  idx = (traj-traj_start)/traj_incr/bin_size
  data_binned[:,0,idx] = data[:,0,idx*bin_size]
  data_binned[:,1,idx] = data[:,1,idx*bin_size]
  data_binned[:,2:,idx] = np.average(data[:,2:,idx*bin_size:(idx+1)*bin_size], axis=2)
  fout = data_stem + ".bin_" + str(bin_size) + "." + str(traj)
  f = open(fout,'w')
  for row in xrange(0,Nrows):
    line = "{0:3d}\t{1:3d}".format(int(data_binned[row,0,idx]),int(data_binned[row,1,idx]))
    for col in xrange(2,Ncols):
      line += "\t{0:1.10e}".format(data_binned[row,col,idx])
    print(line, file=f)
  f.close()

print("done.")
