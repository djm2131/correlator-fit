#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "fitter.h"
#include "fitter_controls.h"
#include "xml_parser.h"
#include "jack_stats.h"

// Data to pass to GSL routines through void*
typedef struct {
  double* t;
  double* C;
  Fitter* me;
} fit_data;

Fitter::Fitter(std::string xml_path)
{
  printf("\n----- Constructing fitter and loading raw data -----\n");
  
  {
    XML_parser XML(xml_path);
    fc = XML.parse_all();
  }
  
  corrs.resize(fc.fits.size());
  for(unsigned int i=0; i<fc.fits.size(); ++i){
    corrs[i] = new Correlator(fc, i);
  }
  
  printf("\n----- done -----\n");
}

int Fitter::f(const gsl_vector* x, void* data, gsl_vector* y)
{
  double* t = ((fit_data*)data)->t;
  double* C = ((fit_data*)data)->C;
  
  int p_idx(0), y_idx(0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    // Get the current parameter values for this correlator
    int Np = corrs[i]->get_Np();
    std::vector<double> p(Np);
    for(int j=0; j<Np; ++j){ p[j] = gsl_vector_get(x, p_idx); ++p_idx; }
    
    // Compute difference between fit and data for each data pt.
    int Ndat = fc.fits[i].t_max - fc.fits[i].t_min + 1;
    for(int j=0; j<Ndat; ++j){
      gsl_vector_set(y, y_idx, corrs[i]->eval(t[y_idx],p) - C[y_idx]);
      ++y_idx;
    }
  }
  
  return GSL_SUCCESS;
}

// Neat trick: we can construct a C function pointer to a member
// function of a C++ class by declaring a static wrapper which
// gets passed a 'this' pointer through the fit_data struct
int Fitter::f_wrapper(const gsl_vector* x, void* data, gsl_vector* y)
{
  return static_cast<fit_data*>((fit_data*)data)->me->f(x,data,y);
}

int Fitter::df(const gsl_vector* x, void* data, gsl_matrix* J)
{
  double *t = ((fit_data*)data)->t;
  
  int p_idx(0), y_idx(0), k_offset(0);
  for(unsigned int i=0; i<corrs.size(); ++i)
  {
    int Np = corrs[i]->get_Np();
    std::vector<double> p(Np), df(Np);
    for(int j=0; j<Np; ++j){ p[j] = gsl_vector_get(x, p_idx); ++p_idx; }
    
    int Ndat = fc.fits[i].t_max - fc.fits[i].t_min + 1;
    for(int j=0; j<Ndat; ++j){
      df = corrs[i]->eval_derivs(t[y_idx],p);
      for(int k=0; k<fc.Nparams; ++k){ 
        if((k >= k_offset) && (k <= k_offset+Np)){ gsl_matrix_set(J, y_idx, k, df[k-k_offset]); }
        else{ gsl_matrix_set(J, y_idx, k, 0.0); }
      }
      ++y_idx;
    }
  }
  
  return GSL_SUCCESS;
}

int Fitter::df_wrapper(const gsl_vector* x, void* data, gsl_matrix* J)
{
  return static_cast<fit_data*>((fit_data*)data)->me->df(x,data,J);
}

fit_results Fitter::do_fit(void)
{
  std::vector<double> chi2pdof_jacks(fc.Ntraj);
  fit_results fr = { 0.0, 0.0, std::vector<double>(fc.Nparams), std::vector<double>(fc.Nparams), 
                     std::vector<std::vector<double>>(fc.Ntraj,std::vector<double>(fc.Nparams)) };
  
  // Loop over jackknife samples
  // -1 is the fit to all data (central value)
  for(int jknife_idx=-1; jknife_idx<fc.Ntraj; ++jknife_idx)
  {
    if(jknife_idx == -1){ printf("\n----- Fitting to all data -----\n"); }
    else{ printf("\n----- Jackknife sample %d of %d -----\n", jknife_idx+1, fc.Ntraj); }
    
    double* weights = new double[fc.Ndata];
    fit_data fd;
    fd.t = new double[fc.Ndata];
    fd.C = new double[fc.Ndata];
    fd.me = this;
    {
      int ii = 0;
      std::vector<double> Cavg(fc.Ndata), Cstd(fc.Ndata);
      std::vector<std::vector<double>> Cfit(fc.Ntraj, std::vector<double>(fc.Ndata,0.0));
      
      // Loop through the raw data and store the points to fit
      for(unsigned int i=0; i<fc.fits.size(); ++i)
      {
        // t
        int t_min = fc.fits[i].t_min;
        int t_max = fc.fits[i].t_max;
        int t_len = t_max - t_min + 1;
        for(int j=0; j<t_len; ++j){ 
          fd.t[ii] = t_min+j; 
          ii += 1; 
        }
        
        // C(t)
        for(int j=0; j<fc.Ntraj; ++j){
          int l = 0;
          for(int k=0; k<corrs[i]->get_corr_ndata(); ++k){
            if(corrs[i]->include_data_pt(j,k)){ 
              Cfit[j][l] = corrs[i]->get_data_pt(j,k);
              l += 1;
            }
          }
        }
      }
      
      if(jknife_idx >= 0){ Cfit.erase(Cfit.begin()+jknife_idx); }
      
      // Compute jackknife averages and errors
      Cavg = jack_avg(Cfit);
      Cstd = jack_std(Cfit, Cavg, true);
      for(int i=0; i<fc.Ndata; ++i){ 
        fd.C[i] = Cavg[i]; 
        weights[i] = pow(Cstd[i],-2.0); 
      }
    }
    
    // Initialize GSL objects
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    gsl_matrix *J = gsl_matrix_alloc(fc.Ndata, fc.Nparams);
    gsl_multifit_function_fdf f;
    gsl_vector *p = gsl_vector_alloc(fc.Nparams);
    gsl_vector_view w = gsl_vector_view_array(weights, fc.Ndata);
    gsl_vector *res_f;
    double chi, chi0;
    
    // Set initial params to user-supplied guesses
    int p_idx(0);
    for(unsigned int i=0; i<corrs.size(); ++i){
      for(int j=0; j<corrs[i]->get_Np(); ++j){ 
        gsl_vector_set(p, p_idx, fc.p0[i][j]);
        ++p_idx;
      }
    }

    f.f = &Fitter::f_wrapper;
    (fc.numerical_derivs) ? (f.df = NULL) : (f.df = &Fitter::df_wrapper);
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
    
    // Catch a weird GSL bug --- why does this happen !?
    if((gsl_multifit_fdfsolver_niter(s) == 1) && (info != 1)){ 
      printf("\nError: GSL is being tempermental. Try running the program again.\n\n"); 
      exit(-1); 
    }

    // Summary of fit results
    printf("\n** Summary from method '%s' **\n", gsl_multifit_fdfsolver_name(s));
    printf("Status: %s\n", gsl_strerror(status));
    printf("Number of iterations: %zu\n", gsl_multifit_fdfsolver_niter(s));
    printf("Function evaluations: %zu\n", f.nevalf);
    printf("Jacobian evaluations: %zu\n", f.nevaldf);
    printf("Reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    printf("Initial chisq = %g\n", chi0);
    printf("Final chisq = %g\n", chi);
    { 
      double dof = fc.Ndata - fc.Nparams;
      printf("chisq/dof = %g\n",  pow(chi,2.0)/dof);
      printf("Fit parameters:\n");
      for(int i=0; i<fc.Nparams; ++i){
        (jknife_idx == -1) ? (fr.p_cv[i] = gsl_vector_get(s->x,i)) : (fr.p_jacks[jknife_idx][i] = gsl_vector_get(s->x,i));
        printf("  %s = %1.8e\n", fc.p_names[i].c_str(), gsl_vector_get(s->x,i));
      }
      (jknife_idx == -1) ? (fr.chi2pdof = pow(chi,2.0)/dof) : (chi2pdof_jacks[jknife_idx] = pow(chi,2.0)/dof);
    }

    // Clean up
    gsl_multifit_fdfsolver_free(s);
    gsl_matrix_free(J);
    gsl_vector_free(p);
    
    delete[] fd.C;
    delete[] fd.t;
    delete[] weights;
  }
    
  fr.p_err = jack_std(fr.p_jacks, fr.p_cv, false);
  fr.chi2pdof_err = jack_std(chi2pdof_jacks, fr.chi2pdof, false);
  printf("\n----- done -----\n");
  
  return fr;
}

void Fitter::print_results(fit_results& fr)
{
  printf("\nFit results:\n");
  for(unsigned int i=0; i<fr.p_cv.size(); ++i){
    printf("  %s = %1.8e +/- %1.8e\n", fc.p_names[i].c_str(), fr.p_cv[i], fr.p_err[i]);
  }
  printf("  chi^2/dof = %1.8e +/- %1.8e\n", fr.chi2pdof, fr.chi2pdof_err);
}

Fitter::~Fitter(){ for(unsigned int i=0; i<corrs.size(); ++i){ delete corrs[i]; } }