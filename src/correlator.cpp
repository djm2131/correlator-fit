/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: src/correlator.cpp

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "correlator.h"
#include "fitter_controls.h"

std::vector<double> tokenize_line(const std::string& s, const std::string& delim)
{
  std::istringstream iss(s);
  std::string token;
  std::vector<double> tokens;

  while(getline(iss,token,' ')){ if(!token.empty()){ tokens.push_back( strtod(token.c_str(), nullptr) ); } }
  // if(tokens.size() != 3){
  //   printf("Error: expected data format is Nx3 matrix with format: (t, real(C(t)), imag(C(t))).\n");
  //   exit(-1);
  // }

  return tokens;
}

Correlator::Correlator(fitter_controls& fc, const int& this_fit)
{
  // Add number of data points in this fit to total
  fc.Ndata += static_cast<int>( fc.fits[this_fit].t_max - fc.fits[this_fit].t_min + 1 );

  ////////////////////////////////////////////////////////////////
  // Read correlator data from disk.
  //
  // This assumes a file containing an Nx3 matrix formatted as
  //    t   real(C(t))   imag(C(t))
  // for each trajectory included in the fit.
  // 3-pt functions should be treated as a simultaneous fit of a
  // series of 2-pt correlators for each source-sink separation.
  ////////////////////////////////////////////////////////////////
  int traj = fc.traj_start;
  t.resize(fc.Ntraj);
  C.resize(fc.Ntraj);
  if(!fc.fits[this_fit].resample){ w.resize(fc.Ntraj); }
  include_in_fit.resize(fc.Ntraj);

  for(int i=0; i<fc.Ntraj; ++i)
  {
    // Path to data for this trajectory
    std::stringstream ss;
    ss << fc.fits[this_fit].data_stem << "." << traj;

    // Open data file
    std::ifstream in(ss.str());
    if(!in.is_open()){
      printf("Error: failed to open file %s\n", ss.str().c_str());
      exit(-1);
    }

    // Parse line by line and store result
    std::string line;
    std::vector<double> tokens;
    while(getline(in,line)){
      tokens = tokenize_line(line, " ");
      t[i].push_back(tokens[0]);
      C[i].push_back(tokens[1]);
      if(!fc.fits[this_fit].resample){ w[i].push_back( pow(tokens[2],-2.0) ); }
      if((tokens[0] >= fc.fits[this_fit].t_min) && (tokens[0] <= fc.fits[this_fit].t_max)){
        include_in_fit[i].push_back(1);
      } else {
        include_in_fit[i].push_back(0);
      }
    }

    // some debug messages
#if(0)
    printf("\nConstructing correlator for fit %d:\n", this_fit);
      printf("\n--- Fitting to t = (%d,%d)\n", static_cast<int>(fc.fits[this_fit].t_min),
        static_cast<int>(fc.fits[this_fit].t_max));
      printf("\n--- Got data for trajectory %d:\n", traj);
      for(unsigned int j=0; j<t[i].size(); ++j){
        printf("%3d\t%1.6e\n", static_cast<int>(t[i][j]), C[i][j]);
      }
#endif

    traj += fc.traj_inc;
  }

  // If we aren't resampling, append the results
  // computed using all data to the beginning
  if(!fc.fits[this_fit].resample){
    std::stringstream ss;
    ss << fc.fits[this_fit].data_stem;
    std::ifstream in(ss.str());
    if(!in.is_open()){
      printf("Error: failed to open file %s\n", ss.str().c_str());
      exit(-1);
    }
    std::string line;
    std::vector<double> tokens;
    std::vector<double> t_tmp;
    std::vector<double> C_tmp;
    std::vector<double> w_tmp;
    std::vector<int> include_tmp;
    while(getline(in,line)){
      tokens = tokenize_line(line, " ");
      t_tmp.push_back(tokens[0]);
      C_tmp.push_back(tokens[1]);
      w_tmp.push_back(pow(tokens[2],-2.0));
      if((tokens[0] >= fc.fits[this_fit].t_min) && (tokens[0] <= fc.fits[this_fit].t_max)){
        include_tmp.push_back(1);
      } else {
        include_tmp.push_back(0);
      }
    }
    t.insert(t.begin(), t_tmp);
    C.insert(C.begin(), C_tmp);
    w.insert(w.begin(), w_tmp);
    include_in_fit.insert(include_in_fit.begin(), include_tmp);
  }

  // Get the covariance matrix if we're doing correlated fits
  if(fc.correlated_fits)
  {
    size_t N = static_cast<size_t >( fc.fits[this_fit].t_max - fc.fits[this_fit].t_min + 1 );
    mcov.resize(N, std::vector<double>(N));

    std::stringstream ss;
    ss << fc.fits[this_fit].cov_matrix_stem;
    std::ifstream in(ss.str());
    if(!in.is_open()){
      printf("Error: failed to open file %s\n", ss.str().c_str());
      exit(-1);
    }

    int t1(0), row(0);
    std::string line;
    std::vector<double> tokens;
    while(getline(in,line))
    {
      if((t1 < fc.fits[this_fit].t_min) || (t1 > fc.fits[this_fit].t_max)){ t1++; continue; }
      tokens = tokenize_line(line, " ");
      int col(0);
      for(int t2=0; t2<tokens.size(); t2++){
        if((t2 >= fc.fits[this_fit].t_min) && (t2 <= fc.fits[this_fit].t_max)){
          mcov[row][col] = tokens[t2];
          col++;
        }
      }
      row++;
      t1++;
    }

    gsl_pinv(mcov, fc.svd_cut);
  }

  // Construct the appropriate FitFunc object for this fit form
  // The available forms are defined in include/fit_functions.h
  // double V = static_cast<double>( pow(fc.L,3)*fc.T );
  double V = pow(fc.L, 3);
  if     (fc.fits[this_fit].fit_type == "constant"       ){ f = new ConstFunc();                              }
  else if(fc.fits[this_fit].fit_type == "linear"         ){ f = new LinearFunc();                             }
  else if(fc.fits[this_fit].fit_type == "exp"            ){ f = new ExpFunc(V);                               }
  else if(fc.fits[this_fit].fit_type == "exp_plus_const" ){ f = new ExpPlusConstFunc(V);                      }
  else if(fc.fits[this_fit].fit_type == "double_exp"     ){ f = new DoubleExpFunc(V);                         }
  else if(fc.fits[this_fit].fit_type == "cosh_pp"        ){ f = new CoshFuncPP(fc.T, V);                      }
  else if(fc.fits[this_fit].fit_type == "cosh_ap"        ){ f = new CoshFuncAP(fc.T, V);                      }
  else if(fc.fits[this_fit].fit_type == "cosh_decay"     ){ f = new CoshFuncDecay(fc.T, V);                   }
  else if(fc.fits[this_fit].fit_type == "cosh_two_pion"  ){ f = new CoshFuncTwoPion(fc.T,V);                  }
  else if(fc.fits[this_fit].fit_type == "ratio_two_pion" ){ f = new TwoPionRatioFunc(fc.T,V);                 }
  else if(fc.fits[this_fit].fit_type == "zv"             ){ f = new ZVFunc(fc.T, V, fc.fits[this_fit].t_sep); }
  else if(fc.fits[this_fit].fit_type == "bk"             ){ f = new BKFunc(fc.T, V, fc.fits[this_fit].t_sep); }
  else{
    printf("Error: unrecognized fit type %s\n", fc.fits[this_fit].fit_type.c_str());
    exit(-1);
  }

  // If we are computing the effective mass for this correlator,
  // construct the appropriate EffMassFunc object as well.
  if(fc.fits[this_fit].do_eff_mass){
    if(fc.fits[this_fit].eff_mass_type == "constant"            ){ ef = new ConstEffMassFunc(); }
    else if(fc.fits[this_fit].eff_mass_type == "log"            ){ ef = new LogEffMassFunc(); }
    else if(fc.fits[this_fit].eff_mass_type == "cosh"           ){ ef = new CoshEffMassFunc(); }
    else if(fc.fits[this_fit].eff_mass_type == "sinh" ||
            fc.fits[this_fit].eff_mass_type == "subtracted_sinh"){ ef = new SinhEffMassFunc(); }
    else {
      printf("Error: unrecognized effective mass type %s\n", fc.fits[this_fit].eff_mass_type.c_str());
      exit(-1);
    }
  }
}

Correlator::~Correlator(){ delete f; }
