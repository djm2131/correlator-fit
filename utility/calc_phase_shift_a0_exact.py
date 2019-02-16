#!/usr/bin/python

import sys
import numpy as np
from scipy.integrate import quad
from scipy.misc import factorial
from itertools import product

# Computes the integrand
#   I(t) = e^{q^{2}t} (\frac{\pi}{t})^{3/2} \sum\limits_{\vec{n} \in \mathbb{Z}^{3}} e^{-\vec{n}^{2} \pi^{2}/t}
def integrand(t,q,ncut):
  I = 0.0;
  for nx,ny,nz in product(xrange(-ncut,ncut+1), repeat=3):
    if nx == ny and nx == nz:
      continue
    I += np.exp(-np.pi**2*(nx**2 + ny**2 + nz**2)/t)
  I *= np.exp(q**2*t) * (np.pi/t)**(1.5)
  return I

def calc_ps(q):

  # Order to truncate sums
  ncut = 4
  lcut = 10

  Z00 = 0.0;

  # First term:
  #   \sum\limits_{\vec{n} \in \mathbb{Z}^{3}} e^{-(\vec{n}^{2} - q^{2})}/(\vec{n}^{2} - q^{2})
  for nx,ny,nz in product(xrange(-ncut,ncut+1), repeat=3):
    Z00 += np.exp(-(nx**2 + ny**2 + nz**2 - q**2))/(nx**2 + ny**2 + nz**2 - q**2)

  # Second term:
  #   \sum\limits_{l=0}^{\infty} \frac{\pi^{3/2}}{l-\frac{1}{2}} \frac{q^{2l}}{l!}
  for l in np.arange(0,lcut+1):
    Z00 += np.pi**(1.5) * q**(2*l) / (l-0.5) / factorial(l)

  # Third term:
  #   \int\limits_{0}^{1} dt I(t)
  # where I(t) is computed by the routine "integrand".
  Z00 += quad(integrand, 0.0, 1.0, args=(q,ncut))[0]

  # return 180.0 * (np.arctan(-np.pi**(1.5)*q*np.sqrt(4.0*np.pi)/Z00) + np.rint(q**2)*np.pi) / np.pi  # deg
  return np.arctan(-np.pi**(1.5)*q*np.sqrt(4.0*np.pi)/Z00) + np.rint(q**2)*np.pi                    # rad

if len(sys.argv) != 7:
  print "Usage: python calc_phase_shift_a0_exact.py <L> <mpi> <mpi_jacks> <Epipi> <Epipi_jacks> <out_dir>"
  exit(0)

L = float(sys.argv[1])
mpi = np.genfromtxt(sys.argv[2])
mpi_jack = np.genfromtxt(sys.argv[3])
if np.array_equal(mpi, mpi_jack):
  mpi = mpi[0]
  mpi_jack = mpi_jack[1:]
dEpipi = np.genfromtxt(sys.argv[4])
dEpipi_jack = np.genfromtxt(sys.argv[5])
if np.array_equal(dEpipi, dEpipi_jack):
  dEpipi = dEpipi[0]
  dEpipi_jack = dEpipi_jack[1:]
out_dir = sys.argv[6]

Epipi = dEpipi + 2.0*mpi
Epipi_jack = dEpipi_jack + 2.0*mpi_jack
p2 = 0.25*Epipi**2 - mpi**2
p = np.sqrt(0.25*Epipi**2 - mpi**2)
delta = -calc_ps(0.5*p*L/np.pi)
a0 = mpi*delta/p

njack = len(mpi_jack)
p_jack = np.zeros((njack))
p2_jack = np.zeros((njack))
delta_jack = np.zeros((njack))
a0_jack = np.zeros((njack))
for jack in np.arange(0,njack):
  p2_jack[jack] = 0.25*Epipi_jack[jack]**2 - mpi_jack[jack]**2
  p_jack[jack] = np.sqrt(p2_jack[jack])
  delta_jack[jack] = -calc_ps(0.5*np.sqrt(p2_jack[jack])*L/np.pi)
  a0_jack[jack] = mpi_jack[jack] * delta_jack[jack] / np.sqrt(p2_jack[jack])

print "mpi = {0} +/- {1}".format(mpi, np.sqrt(njack-1)*np.std(mpi_jack))
print "Epipi = {0} +/- {1}".format(Epipi, np.sqrt(njack-1)*np.std(Epipi_jack))
print "a*p = {0} +/- {1}".format(p, np.sqrt(njack-1)*np.std(p_jack))
print "(a*p)^2 = {0} +/- {1}".format(p2, np.sqrt(njack-1)*np.std(p2_jack))
print "delta(p) = {0} +/- {1}".format(delta*180.0/np.pi, np.sqrt(njack-1)*np.std(delta_jack)*180.0/np.pi)
print "mpi*a0 = {0} +/- {1}\n".format(a0, np.sqrt(njack-1)*np.std(a0_jack))

np.savetxt(out_dir + "/p.dat", [p])
np.savetxt(out_dir + "/p_jacks.dat", p_jack)
np.savetxt(out_dir + "/d02.dat", [delta*180.0/np.pi])
np.savetxt(out_dir + "/d02_jacks.dat", delta_jack*180.0/np.pi)
