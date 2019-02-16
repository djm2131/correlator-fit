#!/usr/bin/python3

import sys
import os.path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

###################################################################
## Plot midpoint of ZV ratio as a function of pion separation
##
## David Murphy (djm2131@columbia.edu)
## 12/30/2016
###################################################################

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
  print("Usage: ./plot_zv_midpt_vs_sep.py <options> <data_stem> <path-to-out-file> <T>")
  print("Options:")
  print("  -legendloc <legend location>")
  print("  -ylabel <label for y-axis>")
  print("  -datalabel <legend lable for data pts.>")
  print("  -xlim <x_min> <x_max>")
  print("  -ylim <y_min> <y_max>")
  print("  -tfit <t_min> <t_max>")
  print("  -fit_data <cv> <jacks>")
  print("Note: if -tfit is set the data pts. included in the fit are marked in red.")
  print("Note: if -fit_data is also set the fitted mass and error is added to the plot.")
  exit(0)

argc = len(sys.argv)
if(argc < 4):
  usage()

# Default values
legendloc = 'lower right'
ylabel = ''
data_label = ''
xlim = [ np.nan, np.nan ]
ylim = [ np.nan, np.nan ]
t_fit = [ np.nan, np.nan ]
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
    data_label = sys.argv[arg+1]
    arg += 2
  elif(sys.argv[arg] == '-xlim'):
    xlim[0] = float(sys.argv[arg+1])
    xlim[1] = float(sys.argv[arg+2])
    arg += 3
  elif(sys.argv[arg] == '-ylim'):
    ylim[0] = float(sys.argv[arg+1])
    ylim[1] = float(sys.argv[arg+2])
    arg += 3
  elif(sys.argv[arg] == '-tfit'):
    t_fit[0] = float(sys.argv[arg+1])
    t_fit[1] = float(sys.argv[arg+2])
    arg += 3
  elif(sys.argv[arg] == '-fit_data'):
    fit_data[0] = sys.argv[arg+1]
    fit_data[1] = sys.argv[arg+2]
    arg += 3
  else:
    print("Unrecognized argument: %s", sys.argv[arg])
    usage()

# Required command line inputs
data_stem = sys.argv[argc-3]
fout      = sys.argv[argc-2]
T         = int(sys.argv[argc-1])

ZV = np.zeros((T,3))
for sep in range(0,T):
  ZV[sep,0] = sep
  fs = data_stem + "." + str(sep)
  if os.path.isfile(fs):
    tmp = np.genfromtxt(fs)
    if sep%2 == 0:
      idx = int(sep/2)
      ZV[sep,1] = tmp[idx,1]
      ZV[sep,2] = tmp[idx,2]
    else:
      idx = int((sep-1)/2)
      ZV[sep,1] = 0.5*( tmp[idx,1] + tmp[idx+1,1] )
      ZV[sep,2] = 0.5*( tmp[idx,2] + tmp[idx+1,2] )


# Plot data
if not np.isnan(t_fit[0]):
  idx_ti = np.where(ZV[:,0] == t_fit[0])[0][0]
  if ZV[-1,0] >= t_fit[1]:
    idx_tf = np.where(ZV[:,0] == t_fit[1])[0][0] + 1
  elif ZV[-1,0] == t_fit[1]-1:
    t_fit[1] -= 1
  else:
    print("Dimensions of data do not match -tfit option.")
    exit(-1)
  if data_label:
    plt.errorbar(ZV[0:idx_ti,0], ZV[0:idx_ti,1], yerr=ZV[0:idx_ti,2], marker='o', color='k', ls='none', label=data_label)
  else:
    plt.errorbar(ZV[0:idx_ti,0], ZV[0:idx_ti,1], yerr=ZV[0:idx_ti,2], marker='o', color='k', ls='none')
  plt.errorbar(ZV[idx_ti:idx_tf,0], ZV[idx_ti:idx_tf,1], yerr=ZV[idx_ti:idx_tf,2], marker='o', color='r', ls='none')
  plt.errorbar(ZV[idx_tf:,0], ZV[idx_tf:,1], yerr=ZV[idx_tf:,2], marker='o', color='k', ls='none')
else:
  if data_label:
    plt.errorbar(ZV[:,0], ZV[:,1], yerr=ZV[:,2], marker='o', color='k', ls='none', label=data_label)
  else:
    plt.errorbar(ZV[:,0], ZV[:,1], yerr=ZV[:,2], marker='o', color='k', ls='none')

# Plot fit
if fit_data[0]:
  if np.isnan(t_fit[0]):
    print("error: must supply fit range to use -fit_data option.")
    exit(-1)
  cv = np.genfromtxt(fit_data[0])
  jacks = np.genfromtxt(fit_data[1])
  err = np.sqrt(len(jacks)-1.0)*np.std(jacks, ddof=0)
  t = ZV[idx_ti:idx_tf,0]
  plt.fill_between(t, np.tile(cv-err,(len(t),)), np.tile(cv+err,(len(t),)), color='r', alpha=0.5, label='Fit', zorder=0)

# Limits
if not np.isnan(xlim[0]):
  plt.xlim(xlim[0], xlim[1])
if not np.isnan(ylim[0]):
  plt.ylim(ylim[0], ylim[1])

# Labels and legend
plt.xlabel('$\\left| t_{\\rm src} - t_{\\rm snk} \\right| / a$', fontsize=18)
if ylabel:
  plt.ylabel(ylabel, fontsize=18)
if data_label or fit_data[0]:
  handles, labels = plt.gca().get_legend_handles_labels()
  plt.gca().legend(reversed(handles), reversed(labels), loc=legendloc, numpoints=1)

# Save
plt.savefig(fout, bbox_inches='tight')
