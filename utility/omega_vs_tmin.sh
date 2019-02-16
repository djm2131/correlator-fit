#!/bin/bash

fit_binary=/home/dmurphyiv/Dropbox/LQCD/Projects/correlator-fit/build/correlator-fit
results_dir=/home/dmurphyiv/Dropbox/LQCD/Projects/coarse_ensembles/32ID/ms0p0850/correlator_fits/results/fit_params

for ((tmin_box=3; tmin_box<=6; tmin_box++)); do
  for ((tmin_wall=5; tmin_wall<=8; tmin_wall++)); do
    sed -e s/TMIN_BOX/$tmin_box/g -e s/TMIN_WALL/$tmin_wall/g xml/fit_omega.xml.template > xml/fit_omega.xml
    $fit_binary xml/fit_omega.xml > logs/fit_omega_${tmin_box}_${tmin_wall}.log 2>&1
  done
done
