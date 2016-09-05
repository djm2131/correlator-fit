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
  
  // Lattice volume
  int L;
  int T;
  
  // Trajectories / binning
  int traj_start;
  int traj_end;
  int traj_inc;
  int bin_size;
  
  // Parameters for individual correlator fits
  std::vector<fit_params> fits;
} fitter_controls;

#endif