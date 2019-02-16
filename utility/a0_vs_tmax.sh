#!/bin/bash

L=24
fit_binary=/home/dmurphyiv/Dropbox/LQCD/Projects/correlator-fit/build/correlator-fit
a0_binary=/home/dmurphyiv/Dropbox/LQCD/Projects/correlator-fit/utility/calc_a0_asymptotic_luscher.py
results_dir=/home/dmurphyiv/Dropbox/LQCD/Projects/coarse_ensembles/24ID/correlator_fits/results/fit_params

for ((tmin=6; tmin<=6; tmin++)); do
  for ((tmax=tmin+3; tmax<=20; tmax++)); do
    sed -e s/TMIN/$tmin/g -e s/TMAX/$tmax/g xml/fit_2pion.xml.template > xml/fit_2pion.xml
    $fit_binary xml/fit_2pion.xml > /dev/null 2>&1
    $a0_binary $L $results_dir/mpi.dat $results_dir/mpi_jacks.dat $results_dir/Epipi.dat $results_dir/Epipi_jacks.dat $results_dir > logs/a0_${tmin}_${tmax}.log
  done
done
