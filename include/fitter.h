/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/fitter.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __FITTER_H_INCLUDED__
#define __FITTER_H_INCLUDED__

#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "correlator.h"
#include "fitter_controls.h"
#include "xml_parser.h"
#include "fit_results.h"

class Fitter;

// Data to pass to GSL routines through void*
typedef struct {
  double* t;
  double* C;
  double* w;
  double** mcov;
  Fitter* me;
} fit_data;

class Fitter {
  private:
    fitter_controls fc;
    std::vector<Correlator*> corrs;
    std::vector<std::vector<int>> bind_map;

  public:
    explicit Fitter(const std::string& xml_path);
    int get_pidx(const int& corr_idx, const int& corr_p_idx) const;
    int get_corr_idx(const int& i) const;
    void apply_constraints(const gsl_vector* x, std::vector<double>& p, const int& corr_idx) const;
    bool free_param(const int& i) const;
    int f(const gsl_vector* x, void* data, gsl_vector* f) const;
    double chisq_uncorr(const gsl_vector* x, void* data) const;
    double chisq_corr(const gsl_vector* x, void* data) const;
    double chisq(const gsl_vector* x, void* data) const {
      if(fc.correlated_fits){ return chisq_corr(x, data); }
      else{ return chisq_uncorr(x, data); }
    }
    static int f_wrapper(const gsl_vector* x, void* data, gsl_vector* f);
    static double chisq_wrapper(const gsl_vector* x, void* data);
    int df(const gsl_vector* x, void* data, gsl_matrix* J) const;
    static int df_wrapper(const gsl_vector* x, void* data, gsl_matrix* J);
    double LM_fit(const int& jknife_idx, fit_data& fd, fit_results& fr) const;
    double NM_fit(const int& jknife_idx, fit_data& fd, fit_results& fr) const;
    fit_results do_fit();
    void print_results(const fit_results& fr) const;
    void save_jacks(const fit_results& fr) const;
    void save_chi2pdof(const fit_results& fr) const;
    void compute_eff_mass(fit_results& fr) const;
    void save_eff_mass(const fit_results& fr) const;
    ~Fitter();
};

#endif 
