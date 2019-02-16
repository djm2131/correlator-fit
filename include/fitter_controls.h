/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/fitter_controls.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __FITTER_CONTROLS_H_INCLUDED__
#define __FITTER_CONTROLS_H_INCLUDED__

#include <string>
#include <vector>

typedef struct {
  bool resample;
  bool do_eff_mass;
  bool subtract_ts;
  std::string fit_type;
  std::string eff_mass_type;
  std::string data_stem;
  std::string eff_mass_stem;
  double t_min;
  double t_max;
  double t_sep;
} fit_params;

typedef struct {

  // Fitter parameters
  std::string algorithm;
  int numerical_derivs;
  int max_iter;
  double xtol;
  double gtol;
  double ftol;

  // Lattice volume
  int L;
  int T;

  // Trajectories / binning
  int traj_start;
  int traj_end;
  int traj_inc;
  int Ntraj;
  int bin_size;

  // Fit parameters
  int Ndata;                           // # data points
  int Nparams;                         // # fit parameters (ignoring constraints)
  int Nfreeparams;                     // # free parameters
  std::vector<std::string> p_names;    // parameter names
  std::vector<std::vector<double>> p0; // initial guesses

  // Parameters for individual correlator fits
  std::vector<fit_params> fits;

  // Constraints
  bool constrained_fit;
  std::vector<std::vector<int>> p_bindings;

  // Output options
  bool save_jacks;
  std::string jacks_dir;
} fitter_controls;

#endif 
