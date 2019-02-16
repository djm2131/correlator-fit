#!/usr/bin/python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

a0 = np.zeros((1,4))

idx = 0
with open('foo','r') as f:
  for line in f:
    tmin   = float(line.split('_')[1])
    tmax   = float(line.split('_')[2].split('.')[0])
    a0_cv  = float(line.split()[-3])
    a0_err = float(line.split()[-1])
    if idx == 0:
      a0[idx,:] = np.array([[tmin, tmax, a0_cv, a0_err]])
    else:
      a0 = np.append(a0, np.array([[tmin, tmax, a0_cv, a0_err]]), axis=0)
    idx += 1

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

labels = []
for idx in range(0,a0.shape[0]):
  labels.append('[' + str(int(a0[idx,0])) + ',' + str(int(a0[idx,1])) + ']')
 
plt.fill_between([-10,a0.shape[0]+10], [-0.048,-0.048], [-0.038,-0.038], color='k', alpha=0.5, label='Expt.')
plt.errorbar(np.arange(0,a0.shape[0]), a0[:,2], yerr=a0[:,3], marker='o', color='k', ls='none')
plt.xlim(-1,a0.shape[0])
plt.xticks(np.arange(0,a0.shape[0]), labels, rotation='vertical')
plt.xlabel('$[t_{\\rm min}/a,t_{\\rm max}/a]$', fontsize=18)
plt.ylabel('$m_{\\pi} a_{0}^{2}$', fontsize=18)
plt.gca().legend(loc='lower right', numpoints=1)
plt.savefig('a0_vs_tmax.pdf', bbox_inches='tight')
