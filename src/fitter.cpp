/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: src/fitter.cpp

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>

#include "fitter.h"
#include "fitter_controls.h"
#include "xml_parser.h"
#include "jack_stats.h"

int Fitter::get_pidx(const int& corr_idx, const int& corr_p_idx) const
{
  int p_idx(corr_p_idx);
  for(int i=0; i<corr_idx; ++i){ p_idx += corrs[i]->get_Np(); }
  return p_idx;
}

int Fitter::get_corr_idx(const int& i) const
{
  int current_idx(-1);
  for(unsigned int ii=0; ii<corrs.size(); ++ii){
    for(int j=0; j<corrs[i]->get_Np(); ++j){
      current_idx++;
      if(current_idx == i){ return ii; }
    }
  }
  printf("Error: failed to find correlator.\n");
  exit(0);
}

Fitter::Fitter(const std::string& xml_path)
{
  printf("\n----- Constructing fitter and loading raw data -----\n");

  // Read XML file
  {
    XML_parser XML(xml_path);
    fc = XML.parse_all();
  }

  // Construct correlators and compute # fit params
  {
    corrs.resize(fc.fits.size());
    int N_fit_params(0), N_free_fit_params(0);
    for(unsigned int i=0; i<fc.fits.size(); ++i){
      corrs[i] = new Correlator(fc, i);
      N_fit_params += corrs[i]->get_Np();
    }
    fc.Nparams = N_fit_params;
    for(int i=0; i<N_fit_params; ++i){
      if(free_param(i)){ ++N_free_fit_params; }
    }
    fc.Nfreeparams = N_free_fit_params;
  }

  // If some fit parameters are bound to each other,
  // we need to include these constraints in the fits.
  // bind_map[i] contains the indices of all parameters
  // which are bound (i.e. set equal) to parameter i.
  if(fc.constrained_fit){
    int Nc = static_cast<int>(corrs.size());
    bind_map.resize(fc.Nparams);
    for(int i=0; i<Nc; ++i){
      for(int j=0; j<corrs[i]->get_Np(); ++j){
        int this_p = get_pidx(i,j);
        for(unsigned int k=0; k<fc.p_bindings.size(); ++k){
          if((fc.p_bindings[k][0] == i) && (fc.p_bindings[k][1] == j)){
            int bind_to = get_pidx(fc.p_bindings[k][2], fc.p_bindings[k][3]);
            bind_map[this_p].push_back(bind_to);
          }}
      }}
    std::vector<std::vector<int>> bind_map_copy = bind_map;
    for(unsigned int i=0; i<bind_map_copy.size(); ++i){
      for(unsigned int j=0; j<bind_map_copy[i].size(); ++j){
        bind_map[bind_map_copy[i][j]].push_back(i);
        for(unsigned int k=0; k<bind_map_copy[i].size(); ++k){
          if(j == k){ continue; }
          bind_map[bind_map_copy[i][j]].push_back(bind_map_copy[i][k]);
        }
      }}
  }

  /*if(fc.constrained_fit){
    for(int i=0; i<fc.Nparams; ++i){
    if(free_param(i)){
      fc.Nfreeparams += 1;
    }}
  } else {
    fc.Nfreeparams = fc.Nparams;
  }*/

  printf("\n----- done -----\n");
}

void Fitter::apply_constraints(const gsl_vector* x, std::vector<double>& p, const int& corr_idx) const
{
  for(int i=0; i<corrs[corr_idx]->get_Np(); ++i){
    for(unsigned int j=0; j<fc.p_bindings.size(); ++j){
      if((fc.p_bindings[j][2] == corr_idx) && (fc.p_bindings[j][3] == i)){
        p[i] = gsl_vector_get(x, get_pidx(fc.p_bindings[j][0],fc.p_bindings[j][1]));
        break;
      }
    }}
}

bool Fitter::free_param(const int& i) const
{
  for(unsigned int j=0; j<fc.p_bindings.size(); ++j){
    if(i == get_pidx(fc.p_bindings[j][2],fc.p_bindings[j][3])){ return false; }
  }
  return true;
}

int Fitter::f(const gsl_vector* x, void* data, gsl_vector* y) const
{
  double* t = static_cast<fit_data*>(data) -> t;
  double* C = static_cast<fit_data*>(data) -> C;

  int p_idx(0), y_idx(0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    // Get the current parameter values for this correlator
    size_t Np = corrs[i]->get_Np();
    std::vector<double> p(Np);
    for(int j=0; j<Np; ++j){
      p[j] = gsl_vector_get(x, p_idx);
      if(fc.constrained_fit){ apply_constraints(x, p, i); }
      ++p_idx;
    }

    // Compute difference between fit and data for each data pt.
    int Ndat = static_cast<int>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
    for(int j=0; j<Ndat; ++j){
      gsl_vector_set(y, y_idx, corrs[i]->eval(t[y_idx],p) - C[y_idx]);
      ++y_idx;
    }
  }

  return GSL_SUCCESS;
}

double Fitter::chisq_uncorr(const gsl_vector* x, void* data) const
{
  double* t = static_cast<fit_data*>(data) -> t;
  double* C = static_cast<fit_data*>(data) -> C;
  double* w = static_cast<fit_data*>(data) -> w;

  int p_idx(0), y_idx(0);
  double chi2(0.0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    // Get the current parameter values for this correlator
    size_t Np = corrs[i]->get_Np();
    std::vector<double> p(Np);
    for(int j=0; j<Np; ++j){
      p[j] = gsl_vector_get(x, p_idx);
      if(fc.constrained_fit){ apply_constraints(x, p, i); }
      ++p_idx;
    }

    // Sum chi^2
    int Ndat = static_cast<int>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
    for(int j=0; j<Ndat; ++j){
      chi2 += pow( corrs[i]->eval(t[y_idx],p) - C[y_idx], 2.0) * w[y_idx];
      ++y_idx;
    }
  }

  return chi2;
}

double Fitter::chisq_corr(const gsl_vector* x, void* data) const
{
  double* t     = static_cast<fit_data*>(data) -> t;
  double* C     = static_cast<fit_data*>(data) -> C;
  double** mcov = static_cast<fit_data*>(data) -> mcov;

  int p_idx(0), y_idx(0);
  double chi2(0.0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    // Get the current parameter values for this correlator
    size_t Np = corrs[i]->get_Np();
    std::vector<double> p(Np);
    for(int j=0; j<Np; ++j){
      p[j] = gsl_vector_get(x, p_idx);
      if(fc.constrained_fit){ apply_constraints(x, p, i); }
      ++p_idx;
    }

    // Get differences between data and fit
    int Ndat = static_cast<int>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
    std::vector<double> dy(Ndat);
    for(int j=0; j<Ndat; j++){
      dy[j] = corrs[i]->eval(t[y_idx], p) - C[y_idx];
      ++y_idx;
    }

    // Sum chi^2
    for(int j=0; j<Ndat; ++j){
    for(int k=0; k<Ndat; ++k){
      chi2 += dy[j] * mcov[i][j*Ndat+k] * dy[k];
    }}
  }

  return chi2;
}

// Neat trick: we can construct a C function pointer to a member
// function of a C++ class by declaring a static wrapper which
// gets passed a 'this' pointer through the fit_data struct
int Fitter::f_wrapper(const gsl_vector* x, void* data, gsl_vector* y) {
  return static_cast<fit_data*>(data) -> me -> f(x, data, y);
}

double Fitter::chisq_wrapper(const gsl_vector* x, void* data) {
  return static_cast<fit_data*>(data) -> me -> chisq(x, data);
}

int Fitter::df(const gsl_vector* x, void* data, gsl_matrix* J) const
{
  double *t = static_cast<fit_data*>(data) -> t;

  int p_idx(0), y_idx(0), k_offset(0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    size_t Np = corrs[i]->get_Np();
    std::vector<double> p(Np), df(Np);
    for(int j=0; j<Np; ++j){
      p[j] = gsl_vector_get(x, p_idx);
      if(fc.constrained_fit){ apply_constraints(x, p, i); }
      ++p_idx;
    }

    size_t Ndat = static_cast<size_t>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
    for(int j=0; j<Ndat; ++j){
      df = corrs[i]->eval_derivs(t[y_idx],p);
      for(int k=0; k<fc.Nparams; ++k){
        if((k >= k_offset) && (k < k_offset+Np)){ gsl_matrix_set(J, y_idx, k, df[k-k_offset]); }
        else{ gsl_matrix_set(J, y_idx, k, 0.0); }
      }
      ++y_idx;
    }

    // Apply constraints to derivatives:
    // we should have df/dp_i = df/dp_j
    // for any bound parameters p_i and p_j
    if(fc.constrained_fit){
      for(int i=0; i<fc.Ndata; ++i){
        for(int j=0; j<fc.Nparams; ++j){
          if(gsl_matrix_get(J,i,j)!= 0){
            for(unsigned int k=0; k<bind_map[j].size(); ++k){
              gsl_matrix_set(J, i, bind_map[j][k], gsl_matrix_get(J,i,j));
            }}
        }}
    }

    k_offset += Np;
  }

  return GSL_SUCCESS;
}

int Fitter::df_wrapper(const gsl_vector* x, void* data, gsl_matrix* J) {
  return static_cast<fit_data*>(data) -> me -> df(x,data,J);
}

double Fitter::LM_fit(const int& jknife_idx, fit_data& fd, fit_results& fr) const
{
  // Initialize GSL objects
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* s;
  int status, info, dof;
  gsl_matrix *J = gsl_matrix_alloc(fc.Ndata, fc.Nparams);
  gsl_multifit_function_fdf f;
  gsl_vector *p = gsl_vector_alloc(fc.Nparams);
  gsl_vector_view w = gsl_vector_view_array(fd.w, fc.Ndata);
  gsl_vector *res_f;
  double chi(0.0), chi0(0.0);

  // Set initial params to user-supplied guesses
  int p_idx(0);
  for(unsigned int i=0; i<corrs.size(); ++i){
    for(int j=0; j<corrs[i]->get_Np(); ++j){
      gsl_vector_set(p, p_idx, fc.p0[i][j]);
      ++p_idx;
    }
  }

  // Apply constraints to initial guesses
  if(fc.constrained_fit){
    for(unsigned int i=0; i<fc.p_bindings.size(); ++i){
      int pidx_1 = get_pidx(fc.p_bindings[i][0],fc.p_bindings[i][1]);
      int pidx_2 = get_pidx(fc.p_bindings[i][2],fc.p_bindings[i][3]);
      double val = gsl_vector_get(p, pidx_1);
      gsl_vector_set(p, pidx_2, val);
    }
  }

  f.f = &Fitter::f_wrapper;
  (fc.numerical_derivs) ? (f.df = nullptr) : (f.df = &Fitter::df_wrapper);
  f.n = fc.Ndata;
  f.p = fc.Nparams;
  f.params = &fd;

  // Allocate and initialize GSL Levenberg-Marquardt solver
  s = gsl_multifit_fdfsolver_alloc(T, fc.Ndata, fc.Nparams);
  gsl_multifit_fdfsolver_wset(s, &f, p, &w.vector);

  // Compute initial residual norm
  res_f = gsl_multifit_fdfsolver_residual(s);
  chi0 = gsl_blas_dnrm2(res_f);

  // Solve the system
  status = gsl_multifit_fdfsolver_driver(s, fc.max_iter, fc.xtol, fc.gtol, fc.ftol, &info);
  gsl_multifit_fdfsolver_jac(s, J);

  // Compute final residual
  chi = gsl_blas_dnrm2(res_f);

  // Summary of fit results
  #ifdef USE_OMP
  #pragma omp critical
  #endif
  {
    printf("\n** Summary from method '%s' **\n", gsl_multifit_fdfsolver_name(s));
    printf("Status: %s\n", gsl_strerror(status));
    printf("Number of iterations: %zu\n", gsl_multifit_fdfsolver_niter(s));
    printf("Function evaluations: %zu\n", f.nevalf);
    printf("Jacobian evaluations: %zu\n", f.nevaldf);
    printf("Reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    printf("Initial chisq = %g\n", chi0);
    printf("Final chisq = %g\n", chi);
    dof = fc.Ndata - fc.Nfreeparams;
    printf("chisq/dof = %g\n",  pow(chi,2.0)/dof);
    printf("Fit parameters:\n");
    printf("N_fit_parameters = %d\n", fc.Nfreeparams);
    for(int i=0; i<fc.Nparams; ++i){
      if(fc.constrained_fit){
        if(free_param(i)){
          (jknife_idx == -1) ? (fr.p_cv[i] = gsl_vector_get(s->x,i)) : (fr.p_jacks[jknife_idx][i] = gsl_vector_get(s->x,i));
          printf("  %s = %1.8e\n", fc.p_names[i].c_str(), gsl_vector_get(s->x,i));
        }
      } else {
        (jknife_idx == -1) ? (fr.p_cv[i] = gsl_vector_get(s->x,i)) : (fr.p_jacks[jknife_idx][i] = gsl_vector_get(s->x,i));
        printf("  %s = %1.8e\n", fc.p_names[i].c_str(), gsl_vector_get(s->x,i));
      }
    }
    if(status){ printf("\nError: boot %d failed to converge.\n", jknife_idx); exit(-1); }
  }

  // Clean up
  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(J);
  gsl_vector_free(p);

  return pow(chi,2.0) / static_cast<double>(dof);
}

double Fitter::NM_fit(const int& jknife_idx, fit_data& fd, fit_results& fr) const
{
  // Initialize GSL objects
  const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer* s;
  int status, dof;
  double size, chisq0(0.0), chisq(0.0);

  gsl_multimin_function f;
  gsl_vector* p = gsl_vector_alloc(fc.Nparams);
  gsl_vector* ss = gsl_vector_alloc(fc.Nparams);

  // Set initial params to user-supplied guesses
  int p_idx(0);
  for(unsigned int i=0; i<corrs.size(); ++i){
    for(int j=0; j<corrs[i]->get_Np(); ++j){
      gsl_vector_set(p, p_idx, fc.p0[i][j]);
      gsl_vector_set(ss, p_idx, 0.1*fc.p0[i][j]);
      ++p_idx;
    }
  }

  // Apply constraints to initial guesses
  if(fc.constrained_fit){
    for(unsigned int i=0; i<fc.p_bindings.size(); ++i){
      int pidx_1 = get_pidx(fc.p_bindings[i][0],fc.p_bindings[i][1]);
      int pidx_2 = get_pidx(fc.p_bindings[i][2],fc.p_bindings[i][3]);
      double val = gsl_vector_get(p, pidx_1);
      gsl_vector_set(p, pidx_2, val);
    }
  }

  f.f = &Fitter::chisq_wrapper;
  f.n = fc.Nparams;
  f.params = &fd;

  // Allocate and initialize GSL Nelder-Mead minimizer
  s = gsl_multimin_fminimizer_alloc(T, fc.Nparams);
  gsl_multimin_fminimizer_set(s, &f, p, ss);

  // Initial chi^2
  status = gsl_multimin_fminimizer_iterate(s);
  chisq0 = s->fval;

  // Perform the fit
  int iter(1);
  do
  {
    status = gsl_multimin_fminimizer_iterate(s);

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, fc.ftol);
    chisq = s->fval;

    iter++;
  }
  while(status == GSL_CONTINUE && iter < fc.max_iter);

  chisq = s->fval;

  // Summary of fit results
  #ifdef USE_OMP
  #pragma omp critical
  #endif
  {
    printf("\n** Summary from method '%s' **\n", gsl_multimin_fminimizer_name(s));
    printf("Status: %s\n", gsl_strerror(status));
    printf("Number of iterations: %d\n", iter);
    printf("Initial chisq = %g\n", chisq0);
    printf("Final chisq = %g\n", chisq);
    dof = fc.Ndata - fc.Nfreeparams;
    printf("chisq/dof = %g\n",  chisq/dof);
    printf("Fit parameters:\n");
    printf("N_fit_parameters = %d\n", fc.Nfreeparams);
    for(int i=0; i<fc.Nparams; ++i){
      if(fc.constrained_fit){
        if(free_param(i)){
          (jknife_idx == -1) ? (fr.p_cv[i] = gsl_vector_get(s->x,i)) : (fr.p_jacks[jknife_idx][i] = gsl_vector_get(s->x,i));
          printf("  %s = %1.8e\n", fc.p_names[i].c_str(), gsl_vector_get(s->x,i));
        }
      } else {
        (jknife_idx == -1) ? (fr.p_cv[i] = gsl_vector_get(s->x,i)) : (fr.p_jacks[jknife_idx][i] = gsl_vector_get(s->x,i));
        printf("  %s = %1.8e\n", fc.p_names[i].c_str(), gsl_vector_get(s->x,i));
      }
    }
    if(status){ printf("\nError: boot %d failed to converge.\n", jknife_idx); exit(-1); }
  }

  // Clean up
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(ss);
  gsl_vector_free(p);

  return chisq / static_cast<double>(dof);
}

fit_results Fitter::do_fit()
{
  fit_results fr = { 0.0, 0.0, std::vector<double>(fc.Nparams),
      std::vector<double>(fc.Nparams), std::vector<double>(fc.Ntraj),
      std::vector<std::vector<double>>(fc.Ntraj,std::vector<double>(fc.Nparams)) };

  // Loop over jackknife samples
  // -1 is the fit to all data (central value)
  #ifdef USE_OMP
  #pragma omp parallel for
  #endif
  for(int jknife_idx=-1; jknife_idx<fc.Ntraj; ++jknife_idx)
  {
    #ifdef USE_OMP
    #pragma omp critical
    #endif
    {
      if(jknife_idx == -1){ printf("\n----- Fitting to all data -----\n"); }
      else{ printf("\n----- Jackknife sample %d of %d -----\n", jknife_idx+1, fc.Ntraj); }
    }

    fit_data fd;
    fd.t    = new double[fc.Ndata];
    fd.C    = new double[fc.Ndata];
    fd.w    = new double[fc.Ndata];
    if(fc.correlated_fits){
      fd.mcov = new double*[fc.fits.size()];
      for(int i=0; i<fc.fits.size(); i++){
        size_t N = static_cast<size_t>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
        fd.mcov[i] = new double[N*N];
      }
    }
    fd.me = this;
    {
      int ii(0);
      std::vector<double> Cavg(fc.Ndata), Cstd(fc.Ndata);

      // Loop through the raw data and store the points to fit
      int this_corr_start_idx(0);
      for(unsigned int i=0; i<fc.fits.size(); ++i)
      {
        // t
        int t_min = static_cast<int>( fc.fits[i].t_min );
        int t_max = static_cast<int>( fc.fits[i].t_max );
        int t_len = t_max - t_min + 1;
        for(int j=0; j<t_len; ++j){
          fd.t[ii] = t_min+j;
          ii += 1;
        }

        // C(t)
        if(fc.fits[i].resample)
        {
          // Get raw data
          std::vector<std::vector<double>> Cfit(fc.Ntraj, std::vector<double>(t_len,0.0));
          for(int j=0; j<fc.Ntraj; ++j){
            int l(0);
            for(int k=0; k<corrs[i]->get_corr_ndata(); ++k){
              if(corrs[i]->include_data_pt(j,k)){
                Cfit[j][l] = corrs[i]->get_data_pt(j,k);
                ++l;
              }
            }
          }

          // Delete one measurement for each jackknife sample
          if(jknife_idx >= 0){ Cfit.erase(Cfit.begin()+jknife_idx); }

          // Compute jackknife average and error
          std::vector<double> Cavg_tmp = jack_avg(Cfit);
          std::vector<double> Cstd_tmp = jack_std(Cfit, Cavg_tmp, true);
          for(int j=0; j<t_len; ++j){
            fd.C[j + this_corr_start_idx] = Cavg_tmp[j];
            fd.w[j + this_corr_start_idx] = pow(Cstd_tmp[j],-2.0);
          }

          this_corr_start_idx += t_len;
        } else {
          int j(jknife_idx+1), l(0);
          for(int k=0; k<corrs[i]->get_corr_ndata(); ++k){
            if(corrs[i]->include_data_pt(j,k)){
              fd.C[l + this_corr_start_idx] = corrs[i]->get_data_pt(j,k);
              fd.w[l + this_corr_start_idx] = corrs[i]->get_weights(j,k);
              ++l;
            }
          }
          this_corr_start_idx += t_len;
        }

        if(fc.correlated_fits){
          for(int i=0; i<fc.fits.size(); i++)
          {
            size_t N = static_cast<size_t>( fc.fits[i].t_max - fc.fits[i].t_min + 1 );
            for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
              fd.mcov[i][N*j+k] = corrs[i]->get_mcov(j,k);
            }}
          }
        }
      }
    }

    double chi2pdof(0.0);
    if(fc.algorithm == "LM"){ chi2pdof = LM_fit(jknife_idx, fd, fr); }
    else if(fc.algorithm == "NM"){ chi2pdof = NM_fit(jknife_idx, fd, fr); }
    else{ printf("Error: unrecognized algorithm %s\n", fc.algorithm.c_str()); exit(-1); }
    (jknife_idx == -1) ? (fr.chi2pdof = chi2pdof) : (fr.chi2pdof_jacks[jknife_idx] = chi2pdof);

    delete[] fd.w;
    delete[] fd.C;
    delete[] fd.t;
    if(fc.correlated_fits){
      for(int i=0; i<fc.fits.size(); i++){ delete[] fd.mcov[i]; }
      delete[] fd.mcov;
    }
  }

  fr.p_err = jack_std(fr.p_jacks, fr.p_cv, false);
  fr.chi2pdof_err = jack_std(fr.chi2pdof_jacks, fr.chi2pdof, false);
  printf("\n----- done -----\n");

  return fr;
}

void Fitter::print_results(const fit_results& fr) const
{
  printf("\nFit results:\n");
  for(int i=0; i<fc.Nparams; ++i){
    if(free_param(i)){
      printf("  %s = %1.8e +/- %1.8e\n", fc.p_names[i].c_str(), fr.p_cv[i], fr.p_err[i]);
    }}
  printf("  chi^2/dof = %1.8e +/- %1.8e\n\n", fr.chi2pdof, fr.chi2pdof_err);
}

void Fitter::save_jacks(const fit_results& fr) const
{
  printf("----- Saving fit parameters to dir %s -----", fc.jacks_dir.c_str());

  for(int i=0; i<fc.Nparams; ++i){
    if(free_param(i)){

      // central value
      std::string fout_cv = fc.jacks_dir + "/" + fc.p_names[i] + ".dat";
      FILE* fcv = fopen(fout_cv.c_str(), "w");
      fprintf(fcv, "%1.8e\n", fr.p_cv[i]);
      fclose(fcv);

      // jackknife samples
      std::string fout_jacks = fc.jacks_dir + "/" + fc.p_names[i] + "_jacks.dat";
      FILE* fj = fopen(fout_jacks.c_str(), "w");
      for(int j=0; j<fc.Ntraj; ++j){
        fprintf(fj, "%1.8e\n", fr.p_jacks[j][i]);
      }
      fclose(fj);

    }}

  printf("\n\n");
}

void Fitter::save_chi2pdof(const fit_results& fr) const
{
  printf("----- Saving chi^2/dof -----");

  // central value
  std::string fout_cv = fc.chi2pdof_stem + ".dat";
  FILE* fcv = fopen(fout_cv.c_str(), "w");
  fprintf(fcv, "%1.8e\n", fr.chi2pdof);
  fclose(fcv);

  // jackknife samples
  std::string fout_jacks = fc.chi2pdof_stem + "_jacks.dat";
  FILE* fj = fopen(fout_jacks.c_str(), "w");
  for(int j=0; j<fc.Ntraj; ++j){
    fprintf(fj, "%1.8e\n", fr.chi2pdof_jacks[j]);
  }
  fclose(fj);

  printf("\n\n");
}

// Computes effective mass and stores this in fit_results structure
void Fitter::compute_eff_mass(fit_results& fr) const
{
  printf("----- Computing effective masses -----\n");

  int Ncorrs = static_cast<int>(corrs.size());
  fr.effm.resize(Ncorrs);

  // Loop over correlators
  int p_offset(0);
  for(int corr_idx=0; corr_idx<Ncorrs; ++corr_idx){
    if(fc.fits[corr_idx].do_eff_mass){
      printf("  type %s for correlator %d...", fc.fits[corr_idx].eff_mass_type.c_str(), corr_idx);

      size_t Ndata = corrs[corr_idx]->get_corr_ndata();
      size_t npts_stencil = corrs[corr_idx]->get_stencil_size();
      size_t npts_eff_mass = Ndata - npts_stencil + 1;
      std::vector<std::vector<double>> eff_mass_jacks(fc.Ntraj+1, std::vector<double>(npts_eff_mass));

      // Initialize effective mass array
      fr.effm[corr_idx].resize(npts_eff_mass);
      for(int j=0; j<npts_eff_mass; ++j){ fr.effm[corr_idx][j].resize(3); }

      // Loop over jackknife samples for this correlator
      for(int jknife_idx=-1; jknife_idx<fc.Ntraj; ++jknife_idx){

        std::vector<double> C_this_jack(Ndata);
        std::vector<std::vector<double>> C(fc.Ntraj+1, std::vector<double>(Ndata));

        // If this is raw data to resample
        if(fc.fits[corr_idx].resample){

          // Get raw data
          std::vector<std::vector<double>> C(fc.Ntraj, std::vector<double>(Ndata));
          for(int j=0; j<fc.Ntraj; ++j){
            for(int k=0; k<Ndata; ++k){
              C[j][k] = corrs[corr_idx]->get_data_pt(j,k);
            }}

          // Delete one measurement for each jackknife sample
          if(jknife_idx >= 0){ C.erase(C.begin()+jknife_idx); }

          // Compute jackknife average
          C_this_jack = jack_avg(C);

        } else { // if data is already resampled

          for(int k=0; k<Ndata; ++k){
            C_this_jack[k] = corrs[corr_idx]->get_data_pt(jknife_idx+1, k);
          }

        }

        // If enabled, subtract out backward propagating state using fit
        if(fc.fits[corr_idx].subtract_ts){
          size_t Np = corrs[corr_idx]->get_Np();
          double t;
          std::vector<double> p(Np);
          for(int k=0; k<Np; ++k){
            int this_p = k + p_offset;
            if(fc.constrained_fit && !bind_map[this_p].empty()){
              this_p = (this_p < bind_map[this_p][0]) ? this_p : bind_map[this_p][0];
            }
            p[k] = (jknife_idx == -1) ? fr.p_cv[this_p] : fr.p_jacks[jknife_idx][this_p];
          }
          for(int j=0; j<Ndata; ++j){
            t = corrs[corr_idx]->get_time_slice(0,j);
            C_this_jack[j] -= corrs[corr_idx]->thermal_state(t,p);
          }
        }

        // Compute effective mass for each stencil
        std::vector<double> stencil(npts_stencil);
        for(int j=0; j<npts_eff_mass; ++j){
          for(int k=0; k<npts_stencil; ++k){ stencil[k] = C_this_jack[j+k]; }
          eff_mass_jacks[jknife_idx+1][j] = corrs[corr_idx]->eff_mass(stencil);
        }
      }

      // Compute jackknifed effective mass
      int start_idx = (npts_stencil%2 == 0) ? (npts_stencil/2-1) : ((npts_stencil-1)/2);
      for(int j=0; j<npts_eff_mass; ++j){
        fr.effm[corr_idx][j][0] = corrs[corr_idx]->get_time_slice(0,j+start_idx);
        fr.effm[corr_idx][j][1] = eff_mass_jacks[0][j];
      }
      eff_mass_jacks.erase(eff_mass_jacks.begin());
      std::vector<double> eff_mass_avg = jack_avg(eff_mass_jacks);
      std::vector<double> eff_mass_err = jack_std(eff_mass_jacks, eff_mass_avg, false);
      for(int j=0; j<npts_eff_mass; ++j){
        fr.effm[corr_idx][j][2] = eff_mass_err[j];
      }

      printf("done.\n");
    }
    p_offset += corrs[corr_idx]->get_Np();
  }
  printf("\n");
}

void Fitter::save_eff_mass(const fit_results& fr) const
{
  printf("----- Saving effective masses -----");

  int Nc = static_cast<int>(corrs.size());
  for(int i=0; i<Nc; ++i){
    if(fc.fits[i].do_eff_mass){
      FILE* f = fopen(fc.fits[i].eff_mass_stem.c_str(), "w");
      for(unsigned int j=0; j<fr.effm[i].size(); ++j){
        fprintf(f, "%3d %1.8e %1.8e\n", static_cast<int>(fr.effm[i][j][0]), fr.effm[i][j][1], fr.effm[i][j][2]);
      }
      fclose(f);
    }}

  printf("\n\n");
}

Fitter::~Fitter(){ for(unsigned int i=0; i<corrs.size(); ++i){ delete corrs[i]; } }
