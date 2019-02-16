#!/usr/bin/python

import sys
import numpy as np
from scipy.optimize import minimize_scalar

def luscher_formula(a0, mpi, Epipi, L):
  return ( Epipi - 2.0*mpi + (4.0*np.pi*a0)/(mpi*L**3)*( 1.0 - 2.837297*(a0/L) + 6.375183*(a0/L)**2 - 8.311951*(a0/L)**3 ) )**2

if len(sys.argv) != 7:
  print "Usage: python calc_a0_asymptotic_luscher.py <L> <mpi> <mpi_jacks> <Epipi> <Epipi_jacks> <out_dir>"

L = float(sys.argv[1])
mpi = np.genfromtxt(sys.argv[2])
mpi_jack = np.genfromtxt(sys.argv[3])
if np.array_equal(mpi, mpi_jack):
  mpi = mpi[0]
  mpi_jack = mpi_jack[1:]
Epipi = np.genfromtxt(sys.argv[4])
Epipi_jack = np.genfromtxt(sys.argv[5])
if np.array_equal(Epipi, Epipi_jack):
  Epipi = Epipi[0]
  Epipi_jack = Epipi_jack[1:]
out_dir = sys.argv[6]

a0 = mpi * minimize_scalar(luscher_formula, bracket=(0.0,2.0), args=(mpi,Epipi,L), tol=1.0e-14).x
a0dL = a0 / mpi / L

njack = len(mpi_jack)
a0_jack = np.zeros((njack))
a0dL_jack = np.zeros((njack))
for jack in np.arange(0,njack):
  a0_jack[jack] = mpi_jack[jack] * minimize_scalar(luscher_formula, bracket=(0.0,2.0), args=(mpi_jack[jack],Epipi_jack[jack],L), tol=1.0e-14).x
  a0dL_jack[jack] = a0_jack[jack] / mpi_jack[jack] / L

print "mpi*a0 = {0} +/- {1}".format(a0, np.sqrt(njack-1)*np.std(a0_jack, axis=0))
print "a0/L = {0} +/- {1}".format(a0dL, np.sqrt(njack-1)*np.std(a0dL_jack, axis=0))

np.savetxt(out_dir + "/a0.dat", [a0])
np.savetxt(out_dir + "/a0_jacks.dat", a0_jack)
