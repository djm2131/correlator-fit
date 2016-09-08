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
  fc = {};
  
  {
    XML_parser XML(xml_path);
    fc = XML.parse_all();
  }
  
  printf("\n----- Constructing fitter and loading raw data -----\n");
  
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
int Fitter::f_wrapper(const gsl_vector* x, void* data, gsl_vector* f)
{
  return static_cast<fit_data*>((fit_data*)data)->me->f(x,data,f);
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

void Fitter::do_fit(void)
{
  // First we fit to the average of all jackknife samples
  double* weights = new double[fc.Ndata];
  fit_data fd;
  fd.t = new double[fc.Ndata];
  fd.C = new double[fc.Ndata];
  fd.me = this;
  {
    int ii = 0;
    std::vector<double> Cavg(fc.Ndata), Cstd(fc.Ndata);
    printf("fc.Ntraj = %d\n", fc.Ntraj);
    printf("fc.Ndata = %d\n", fc.Ndata);
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
    
    // Compute jackknife averages and errors
    Cavg = jack_avg(Cfit);
    Cstd = jack_std(Cfit, Cavg);
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
  gsl_matrix *covar = gsl_matrix_alloc (fc.Nparams, fc.Nparams);
  gsl_multifit_function_fdf f;
  gsl_vector *p = gsl_vector_alloc(fc.Nparams);
  gsl_vector_view w = gsl_vector_view_array(weights, fc.Ndata);
  gsl_vector *res_f;
  double chi, chi0;
  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;
  
  // Set initial params to user-supplied guesses
  int p_idx(0);
  for(unsigned int i=0; i<corrs.size(); ++i){
    for(int j=0; j<corrs[i]->get_Np(); ++j){ 
      gsl_vector_set(p, p_idx, fc.p0[i][j]); 
    }
  }

  f.f = &Fitter::f_wrapper;
  f.df = &Fitter::df_wrapper;
  f.n = fc.Ndata;
  f.p = fc.Nparams;
  f.params = &fd;

  // Allocate and initialize GSL Levenberg-Marquardt solver
  s = gsl_multifit_fdfsolver_alloc(T, fc.Ndata, fc.Nparams);
  gsl_multifit_fdfsolver_wset(s, &f, p, &w.vector);

  // Compute initial residual norm
  res_f = gsl_multifit_fdfsolver_residual(s);
  chi0 = gsl_blas_dnrm2(res_f);

  // Solve the system with a maximum of 1000 iterations
  status = gsl_multifit_fdfsolver_driver(s, 1000, xtol, gtol, ftol, &info);
  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar(J, 0.0, covar);

  // Compute final residual
  chi = gsl_blas_dnrm2(res_f);

  // Summary of fit results
  printf("\nsummary from method '%s'\n", gsl_multifit_fdfsolver_name(s));
  printf("number of iterations: %zu\n", gsl_multifit_fdfsolver_niter(s));
  printf("function evaluations: %zu\n", f.nevalf);
  printf("Jacobian evaluations: %zu\n", f.nevaldf);
  printf("reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
  printf("initial |f(x)| = %g\n", chi0);
  printf("final   |f(x)| = %g\n", chi);
  { 
    double dof = fc.Ndata - fc.Nparams;
    double c = GSL_MAX_DBL(1.0, chi/sqrt(dof)); 
    printf("chisq/dof = %g\n",  pow(chi,2.0)/dof);
    for(int i=0; i<fc.Nparams; ++i){
      printf("p[%d] = %.5f +/- %.5f\n", i, gsl_vector_get(s->x,i), c*sqrt(gsl_matrix_get(covar,i,i)));
    }
  }
  printf("status = %s\n", gsl_strerror(status));

  // Clean up
  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  gsl_matrix_free(J);
  gsl_vector_free(p);
  delete[] fd.C;
  delete[] fd.t;
  delete[] weights;
}

Fitter::~Fitter(){ for(unsigned int i=0; i<corrs.size(); ++i){ delete corrs[i]; } }