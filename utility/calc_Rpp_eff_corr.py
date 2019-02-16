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

def Rpp(App, mp, dEpp, t, T):
  tp = 0.5*T - t
  return App * ( np.cosh(dEpp*tp) + np.sinh(dEpp*tp)/np.tanh(2.0*mp*tp) )

if len(sys.argv) != 9:
  print("Usage: ./calc_Rpp_eff_corr.py <mp_stem> <dEpp_stem> <App_stem> <out_file> <tmin> <tmax> <dt> <T>")
  exit(0)

mp_stem   = sys.argv[1]
dEpp_stem = sys.argv[2]
App_stem  = sys.argv[3]
out_file  = sys.argv[4]
tmin      = float(sys.argv[5])
tmax      = float(sys.argv[6])
dt        = float(sys.argv[7])
T         = float(sys.argv[8])

# Fetch data
mp_cv      = np.genfromtxt( mp_stem   + ".dat"       )
dEpp_cv    = np.genfromtxt( dEpp_stem + ".dat"       )
App_cv     = np.genfromtxt( App_stem  + ".dat"       )
mp_jacks   = np.genfromtxt( mp_stem   + "_jacks.dat" )
dEpp_jacks = np.genfromtxt( dEpp_stem + "_jacks.dat" )
App_jacks  = np.genfromtxt( App_stem  + "_jacks.dat" )

Ntraj      = mp_jacks.shape[0]
Nsamples   = int( (tmax-tmin)/dt + 1 )
Reff_cv    = 0.0
Reff_err   = 0.0
Reff_jacks = np.zeros((Ntraj))
f          = open(out_file, 'w')
for idx in range(0,Nsamples):

  # Compute correlator and error for this t
  t = tmin + idx*dt
  Reff_cv = Rpp(App_cv, mp_cv, dEpp_cv, t, T)
  for jack in range(0,Ntraj):
    App  = App_jacks[jack]
    mp   = mp_jacks[jack]
    dEpp = dEpp_jacks[jack]
    Reff_jacks[jack] = Rpp(App, mp, dEpp, t, T)
  Reff_err = np.sqrt(Ntraj-1.0)*np.std(Reff_jacks, ddof=0)

  # Write to disk
  line = "{0:1.10e} {1:1.10e} {2:1.10e}".format(t, Reff_cv, Reff_err)
  print(line, file=f)

f.close()
