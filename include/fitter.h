#ifndef __FITTER_H_INCLUDED__
#define __FITTER_H_INCLUDED__

#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "correlator.h"
#include "fitter_controls.h"
#include "xml_parser.h"

class Fitter {
private:
  fitter_controls fc;
  std::vector<double> p;
  std::vector<double> p_err;
  std::vector<Correlator*> corrs;
    
public:
  Fitter(std::string xml_path);
  int f(const gsl_vector* x, void* data, gsl_vector* f);
  static int f_wrapper(const gsl_vector* x, void* data, gsl_vector* f);
  int df(const gsl_vector* x, void* data, gsl_matrix* J);
  static int df_wrapper(const gsl_vector* x, void* data, gsl_matrix* J);
  void do_fit(void);
  ~Fitter();
};

#endif