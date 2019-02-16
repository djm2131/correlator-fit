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

def omega_corr(Z, m, Zp, mp, t, V):
  return Z**2/(2.0*m*V)*np.exp(-m*t) + np.sign(Zp)*Zp**2/(2.0*mp*V)*np.exp(-mp*t)

def omega_eff_mass(Z, m, Zp, mp, t, V):
  # return m - np.sign(Zp)*m*(m-mp)*Zp**2 / ( mp*Z**2*np.exp((mp-m)*t) + np.sign(Zp)*m*Zp**2 )
  return np.arcsinh( 0.5*(omega_corr(Z,m,Zp,mp,t-1.0,V)-omega_corr(Z,m,Zp,mp,t+1.0,V))/omega_corr(Z,m,Zp,mp,t,V) )

if len(sys.argv) != 10:
  print("Usage: ./calc_omega_eff_mass.py <Z_stem> <m_stem> <Zp_stem> <mp_stem> <out_file> <tmin> <tmax> <dt> <L>")
  exit(0)

Z_stem     = sys.argv[1]
m_stem     = sys.argv[2]
Zp_stem    = sys.argv[3]
mp_stem    = sys.argv[4]
out_file   = sys.argv[5]
tmin       = float(sys.argv[6])
tmax       = float(sys.argv[7])
dt         = float(sys.argv[8])
L          = float(sys.argv[9])

# Fetch data
Z_cv     = np.genfromtxt( Z_stem  + ".dat"       )
m_cv     = np.genfromtxt( m_stem  + ".dat"       )
Zp_cv    = np.genfromtxt( Zp_stem + ".dat"       )
mp_cv    = np.genfromtxt( mp_stem + ".dat"       )
Z_jacks  = np.genfromtxt( Z_stem  + "_jacks.dat" )
m_jacks  = np.genfromtxt( m_stem  + "_jacks.dat" )
Zp_jacks = np.genfromtxt( Zp_stem + "_jacks.dat" )
mp_jacks = np.genfromtxt( mp_stem + "_jacks.dat" )

Ntraj      = Z_jacks.shape[0]
Nsamples   = int( (tmax-tmin)/dt + 1 )
meff_cv    = 0.0
meff_err   = 0.0
V          = L**3
meff_jacks = np.zeros((Ntraj))
f          = open(out_file, 'w')
for idx in range(0,Nsamples):

  # Compute correlator and error for this t
  t = tmin + idx*dt
  meff_cv = omega_eff_mass(Z_cv, m_cv, Zp_cv, mp_cv, t, V)
  for jack in range(0,Ntraj):
    Z  = Z_jacks[jack]
    m  = m_jacks[jack]
    Zp = Zp_jacks[jack]
    mp = mp_jacks[jack]
    meff_jacks[jack] = omega_eff_mass(Z, m, Zp, mp, t, V)
  meff_err = np.sqrt(Ntraj-1.0)*np.std(meff_jacks, ddof=0)

  # Write to disk
  line = "{0:1.10e} {1:1.10e} {2:1.10e}".format(t, meff_cv, meff_err)
  print(line, file=f)

f.close()
