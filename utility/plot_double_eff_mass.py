#!/usr/bin/python3

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

###################################################################
## Plot effective mass function using outputs from correlator-fit.
## This routine is for the case of a single mass fit to two
## different correlators.
##
## David Murphy (djm2131@columbia.edu)
## 01/06/2017
###################################################################

def find_nearest(array,value):
  return (np.abs(array-value)).argmin()

plt.style.use('seaborn-pastel')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 16

def usage():
  print("Usage: ./plot_eff_mass.py <options> <path-to-data-1> <path-to-data-2> <path-to-out-file>")
  print("Options:")
  print("  -legendloc <legend location>")
  print("  -ylabel <label for y-axis>")
  print("  -datalabel <corr_idx> <legend lable for data pts.>")
  print("  -xlim <x_min> <x_max>")
  print("  -ylim <y_min> <y_max>")
  print("  -tfit <corr_idx> <t_min> <t_max>")
  print("  -fit_data <corr_1> <corr_2>")
  print("Note: if -tfit is set the data pts. included in the fit are marked in red.")
  print("Note: if -fit_data is also set the fitted mass and error is added to the plot.")
  exit(0)

argc = len(sys.argv)
if(argc < 4):
  usage()

# Default values
legendloc = 'lower right'
ylabel = ''
data_label = [ '', '' ]
xlim = [ np.nan, np.nan ]
ylim = [ np.nan, np.nan ]
t_fit = [ [ np.nan, np.nan ], [ np.nan, np.nan ] ]
fit_data = [ '', '' ]

# Parse command line options
arg = 1
while(arg < argc - 3):
  if(sys.argv[arg] == '-legendloc'):
    legendloc = sys.argv[arg+1]
    arg += 2
  elif(sys.argv[arg] == '-ylabel'):
    ylabel = sys.argv[arg+1]
    arg += 2
  elif(sys.argv[arg] == '-datalabel'):
    idx             = int(sys.argv[arg+1])
    data_label[idx] = sys.argv[arg+2]
    arg += 3
  elif(sys.argv[arg] == '-xlim'):
    xlim[0] = float(sys.argv[arg+1])
    xlim[1] = float(sys.argv[arg+2])
    arg += 3
  elif(sys.argv[arg] == '-ylim'):
    ylim[0] = float(sys.argv[arg+1])
    ylim[1] = float(sys.argv[arg+2])
    arg += 3
  elif(sys.argv[arg] == '-tfit'):
    idx           = int(sys.argv[arg+1])
    t_fit[idx][0] = float(sys.argv[arg+2])
    t_fit[idx][1] = float(sys.argv[arg+3])
    arg += 4
  elif(sys.argv[arg] == '-fit_data'):
    fit_data[0] = sys.argv[arg+1]
    fit_data[1] = sys.argv[arg+2]
    arg += 3
  else:
    print("Unrecognized argument: %s", sys.argv[arg])
    usage()

# Required command line inputs
meff_1 = np.genfromtxt(sys.argv[argc-3])
meff_2 = np.genfromtxt(sys.argv[argc-2])
fout = sys.argv[argc-1]

# Plot data (correlator 1)
if not np.isnan(t_fit[0][0]):
  idx_ti = np.where(meff_1[:,0] == t_fit[0][0])[0][0]
  if meff_1[-1,0] >= t_fit[0][1]:
    idx_tf = np.where(meff_1[:,0] == t_fit[0][1])[0][0] + 1
  elif meff_1[-1,0] == t_fit[0][1]-1:
    t_fit[0][1] -= 1
  else:
    print("Dimensions of data do not match -tfit option.")
    exit(-1)
  if data_label[0]:
    plt.errorbar(meff_1[0:idx_ti,0], meff_1[0:idx_ti,1], yerr=meff_1[0:idx_ti,2], marker='o', color='k', ls='none', label=data_label[0])
  else:
    plt.errorbar(meff_1[0:idx_ti,0], meff_1[0:idx_ti,1], yerr=meff_1[0:idx_ti,2], marker='o', color='k', ls='none')
  plt.errorbar(meff_1[idx_ti:idx_tf,0], meff_1[idx_ti:idx_tf,1], yerr=meff_1[idx_ti:idx_tf,2], marker='o', color='r', ls='none')
  plt.errorbar(meff_1[idx_tf:,0], meff_1[idx_tf:,1], yerr=meff_1[idx_tf:,2], marker='o', color='k', ls='none')
else:
  if data_label[0]:
    plt.errorbar(meff_1[:,0], meff_1[:,1], yerr=meff_1[:,2], marker='o', color='k', ls='none', label=data_label[0])
  else:
    plt.errorbar(meff_1[:,0], meff_1[:,1], yerr=meff_1[:,2], marker='o', color='k', ls='none')

# Plot data (correlator 2)
if not np.isnan(t_fit[1][0]):
  idx_ti = np.where(meff_1[:,0] == t_fit[1][0])[0][0]
  if meff_1[-1,0] >= t_fit[1][1]:
    idx_tf = np.where(meff_1[:,0] == t_fit[1][1])[0][0] + 1
  elif meff_1[-1,0] == t_fit[1][1]-1:
    t_fit[1][1] -= 1
  else:
    print("Dimensions of data do not match -tfit option.")
    exit(-1)
  if data_label[1]:
    plt.errorbar(meff_2[0:idx_ti,0], meff_2[0:idx_ti,1], yerr=meff_2[0:idx_ti,2], marker='s', color='b', ls='none', label=data_label[1])
  else:
    plt.errorbar(meff_2[0:idx_ti,0], meff_2[0:idx_ti,1], yerr=meff_2[0:idx_ti,2], marker='s', color='b', ls='none')
  plt.errorbar(meff_2[idx_ti:idx_tf,0], meff_2[idx_ti:idx_tf,1], yerr=meff_2[idx_ti:idx_tf,2], marker='s', color='r', ls='none')
  plt.errorbar(meff_2[idx_tf:,0], meff_2[idx_tf:,1], yerr=meff_2[idx_tf:,2], marker='s', color='b', ls='none')
else:
  if data_label[1]:
    plt.errorbar(meff_2[:,0], meff_2[:,1], yerr=meff_2[:,2], marker='s', color='b', ls='none', label=data_label[1])
  else:
    plt.errorbar(meff_2[:,0], meff_2[:,1], yerr=meff_2[:,2], marker='s', color='b', ls='none')

# Plot fit
if fit_data[0]:
  if np.isnan(t_fit[0][0]):
    print("error: must supply fit range to use -fit_data option.")
    exit(-1)
  meff_fit_1 = np.genfromtxt(fit_data[0])
  meff_fit_2 = np.genfromtxt(fit_data[1])
  t_min = np.min( [ t_fit[0][0], t_fit[1][0] ] )
  t_max = np.max( [ t_fit[0][1], t_fit[1][1] ] )
  t1 = range( find_nearest(meff_fit_1[:,0], t_fit[0][0]), find_nearest(meff_fit_1[:,0], t_fit[0][1])+1 )
  t2 = range( find_nearest(meff_fit_2[:,0], t_fit[1][0]), find_nearest(meff_fit_2[:,0], t_fit[1][1])+1 )
  plt.fill_between(meff_fit_1[t1,0], meff_fit_1[t1,1]-meff_fit_1[t1,2], meff_fit_1[t1,1]+meff_fit_1[t1,2], \
      color='r', alpha=0.5, label='Fit', zorder=0)
  plt.fill_between(meff_fit_2[t2,0], meff_fit_2[t2,1]-meff_fit_2[t2,2], meff_fit_2[t2,1]+meff_fit_2[t2,2], \
      color='r', alpha=0.5, zorder=0)

# Limits
if not np.isnan(xlim[0]):
  plt.xlim(xlim[0], xlim[1])
if not np.isnan(ylim[0]):
  plt.ylim(ylim[0], ylim[1])

# Labels and legend
plt.xlabel('$t/a$', fontsize=18)
if ylabel:
  plt.ylabel(ylabel, fontsize=18)
if data_label or fit_data[0]:
  handles, labels = plt.gca().get_legend_handles_labels()
  plt.gca().legend(reversed(handles), reversed(labels), loc=legendloc, numpoints=1)

# Save
plt.savefig(fout, bbox_inches='tight')
