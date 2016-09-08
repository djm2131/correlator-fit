#ifndef __CORRELATOR_H_INCLUDED__
#define __CORRELATOR_H_INCLUDED__

#include <vector>
#include <string>
#include "fit_functions.h"
#include "fitter_controls.h"

// Parse a line of floating point numbers
// delimited by token into a vector of doubles
std::vector<double> tokenize_line(std::string& s, std::string& token);

class Correlator {
private:
  // std::vector<double> t_fit;
  // std::vector<double> C_fit;
  std::vector<std::vector<int>> include_in_fit;
  std::vector<std::vector<double>> t;
  std::vector<std::vector<double>> C;
  FitFunc *f;

public:  
  Correlator(fitter_controls& fc, int this_fit);
  bool include_data_pt(int i, int j){ return include_in_fit[i][j]; }
  int get_corr_ndata(){ return C[0].size(); } // total # of data pts. for this correlator (includes data not in fit)
  int get_Np(){ return f->get_Np(); }
  double get_data_pt(int i, int j){ return C[i][j]; }
  double eval(double& t, std::vector<double>& p){ return f->eval(t,p); }
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){ return f->eval_derivs(t,p); }
  ~Correlator();
};

#endif