#!/usr/bin/python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mO = np.genfromtxt("mO_vs_tmin.dat")
mOp = np.genfromtxt("mOp_vs_tmin.dat")

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

# Ground state
labels = []
for idx in range(0,mO.shape[0]):
  labels.append('[' + str(int(mO[idx,0])) + ',' + str(int(mO[idx,1])) + ']')
plt.errorbar(np.arange(0,mO.shape[0]), mO[:,2], yerr=mO[:,3], marker='o', color='k', ls='none')
plt.xlim(-1,mO.shape[0])
plt.xticks(np.arange(0,mO.shape[0]), labels, rotation='vertical')
plt.xlabel('$[t_{\\mathrm{min} \\,\\, LZ_{3}B}/a,t_{\\rm{min} \\,\\, LW}/a]$', fontsize=18)
plt.ylabel('$a m_{\\Omega}$', fontsize=18)
plt.savefig('mO_vs_tmin.pdf', bbox_inches='tight')
plt.clf()

# First excited state
labels = []
for idx in range(0,mOp.shape[0]):
  labels.append('[' + str(int(mOp[idx,0])) + ',' + str(int(mOp[idx,1])) + ']')
plt.errorbar(np.arange(0,mOp.shape[0]), mOp[:,2], yerr=mOp[:,3], marker='o', color='k', ls='none')
plt.xlim(-1,mOp.shape[0])
plt.xticks(np.arange(0,mOp.shape[0]), labels, rotation='vertical')
plt.xlabel('$[t_{\\mathrm{min} \\,\\, LZ_{3}B}/a,t_{\\rm{min} \\,\\, LW}/a]$', fontsize=18)
plt.ylabel('$a m_{\\Omega}\'$', fontsize=18)
plt.savefig('mOp_vs_tmin.pdf', bbox_inches='tight')
plt.clf()
