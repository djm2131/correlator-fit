/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)

    Source file: include/eff_mass.h

    Author: David Murphy (djmurphy@mit.edu)

    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __EFF_MASS_H_INCLUDED__
#define __EFF_MASS_H_INCLUDED__

#include <cstdio>
#include <cmath>

class EffMassFunc {

  protected:
    size_t stencil_size;
    std::string EffMassType;

  public:
    EffMassFunc() : stencil_size(0) {};

    const size_t& get_stencil_size() const { return stencil_size; }

    void check_stencil(const std::vector<double>& x) const
    {
      if(x.size() != stencil_size){
        printf("Error: %s expects a %ld-pt stencil.\n", EffMassType.c_str(), stencil_size);
        exit(-1);
      }
    };

    virtual double eff_mass(const std::vector<double>& x) const = 0;

    virtual ~EffMassFunc() = default;
};

class ConstEffMassFunc : public EffMassFunc {
  public:
    ConstEffMassFunc(){ stencil_size = 1; };
    double eff_mass(const std::vector<double>& x) const override { check_stencil(x); return x[0]; };
};

class LogEffMassFunc : public EffMassFunc {
  public:
    LogEffMassFunc(){ stencil_size = 2; };
    double eff_mass(const std::vector<double>& x) const override { check_stencil(x); return log(x[0]/x[1]); };
};

class CoshEffMassFunc : public EffMassFunc {
  public:
    CoshEffMassFunc(){ stencil_size = 3; };
    double eff_mass(const std::vector<double>& x) const override { check_stencil(x); return acosh(0.5*(x[0]+x[2])/x[1]); };
};

class SinhEffMassFunc : public EffMassFunc {
  public:
    SinhEffMassFunc(){ stencil_size = 3; };
    double eff_mass(const std::vector<double>& x) const override { check_stencil(x); return asinh(0.5*(x[0]-x[2])/x[1]); };
};

#endif 
