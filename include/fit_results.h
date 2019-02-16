/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/fit_results.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __FIT_RESULTS_H_INCLUDED__
#define __FIT_RESULTS_H_INCLUDED__

typedef struct {
  double chi2pdof;
  double chi2pdof_err;
  std::vector<double> p_cv;
  std::vector<double> p_err;
  std::vector<std::vector<double>> p_jacks;
  std::vector<std::vector<std::vector<double>>> effm;
} fit_results;

#endif 
