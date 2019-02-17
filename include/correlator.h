/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/correlator.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __CORRELATOR_H_INCLUDED__
#define __CORRELATOR_H_INCLUDED__

#include <vector>
#include <string>

#include "eff_mass.h"
#include "fit_functions.h"
#include "fitter_controls.h"
#include "gsl_pinv.h"

// Parse a line of floating point numbers
// delimited by token into a vector of doubles
std::vector<double> tokenize_line(std::string& s, std::string& token);

class Correlator
{
  private:
    std::vector<std::vector<int>> include_in_fit;
    std::vector<std::vector<double>> t;
    std::vector<std::vector<double>> C;
    std::vector<std::vector<double>> w;
    std::vector<std::vector<double>> mcov;
    FitFunc *f;
    EffMassFunc *ef;

  public:
    Correlator(fitter_controls& fc, const int& this_fit);

    bool include_data_pt(int i, int j){ return static_cast<bool>( include_in_fit[i][j] ); }

    size_t get_corr_njacks() const { return C.size(); }   // # of jackknife samples
    size_t get_corr_ndata() const { return C[0].size(); } // total # of data pts. for this correlator (includes data not in fit)
    size_t get_Np() const { return f->get_Np(); }
    size_t get_stencil_size() const { return ef->get_stencil_size(); }

    const double& get_time_slice(const int& i, const int& j) const { return t[i][j]; }
    const double& get_data_pt(const int& i, const int& j) const { return C[i][j]; }
    const double& get_weights(const int& i, const int& j) const { return w[i][j]; }
    const double& get_mcov(const int& i, const int& j) const { return mcov[i][j]; }

    double eval(double& t, std::vector<double>& p) const { return f->eval(t,p); }
    double thermal_state(double& t, std::vector<double>& p) const { return f->thermal_state(t,p); }
    double eff_mass(std::vector<double>& x) const { return ef->eff_mass(x); }
    std::vector<double> eval_derivs(double& t, std::vector<double>& p) const { return f->eval_derivs(t,p); }

    virtual ~Correlator();
};

#endif 
