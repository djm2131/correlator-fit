#ifndef __FITTER_H_INCLUDED__
#define __FITTER_H_INCLUDED__

#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "correlator.h"
#include "fitter_controls.h"
#include "xml_parser.h"
#include "fit_results.h"

class Fitter {
private:
  fitter_controls fc;
  std::vector<Correlator*> corrs;
  std::vector<std::vector<int>> bind_map;
    
public:
  Fitter(std::string xml_path);
  int get_pidx(int corr_idx, int corr_p_idx);
  void apply_constraints(const gsl_vector* x, std::vector<double>& p, int corr_idx);
  bool free_param(int i);
  int f(const gsl_vector* x, void* data, gsl_vector* f);
  static int f_wrapper(const gsl_vector* x, void* data, gsl_vector* f);
  int df(const gsl_vector* x, void* data, gsl_matrix* J);
  static int df_wrapper(const gsl_vector* x, void* data, gsl_matrix* J);
  fit_results do_fit(void);
  void print_results(fit_results& fr);
  void compute_eff_mass(fit_results& fr);
  void save_eff_mass(fit_results& fr);
  ~Fitter();
};

#endif