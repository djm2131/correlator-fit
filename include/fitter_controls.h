#ifndef __FITTER_CONTROLS_H_INCLUDED__
#define __FITTER_CONTROLS_H_INCLUDED__

#include <string>
#include <vector>

typedef struct {
  std::string fit_type;
  std::string data_stem;
  double t_min;
  double t_max;
} fit_params;

typedef struct {
  
  // Levenberg-Marquardt parameters
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
  int Ndata;
  int Nparams;
  std::vector<std::string> p_names; // parameter names
  std::vector<std::vector<double>> p0; // initial guesses
  
  // Parameters for individual correlator fits
  std::vector<fit_params> fits;
  
  // Constraints
  bool constrained_fit;
  std::vector<std::vector<int>> p_bindings;
} fitter_controls;

#endif