#ifndef __FITFUNC_H_INCLUDED__
#define __FITFUNC_H_INCLUDED__

#include <string>
#include <cmath>

class FitFunc {
    
protected:
  unsigned int Np; // expected number of parameters
  std::string FitType;

public:
  FitFunc(){}
    
  void check_Np(std::vector<double>& p){
    if(p.size() != Np){ Np_error(p.size()); }
  }
    
  void Np_error(unsigned int Ngot){
    if(Np != Ngot){
      printf("Error: %s expects %d parameters but got %d.\n", FitType.c_str(), Np, Ngot);
      exit(-1);
    }
  }
  
  int get_Np(void){ return Np; }
    
  virtual double eval(double& t, std::vector<double>& p) = 0;
    
  virtual std::vector<double> eval_derivs(double& t, std::vector<double>& p) = 0;

  virtual ~FitFunc(){}
};

// Linear function: y[t] = p[0] + p[1]*t
class LinearFunc : public FitFunc{
    
public:    
  LinearFunc() : FitFunc() { Np = 2; FitType = "linear"; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return p[0] + p[1]*t; 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return {1.0, t};
  }
  
  virtual ~LinearFunc(){}
};

// Single exponential: y[t] = p[0]^2/(2.0*p[1]*V) * exp(-p[1]*t)
class ExpFunc : public FitFunc{
  
protected:
  double V;
    
public:    
  ExpFunc(double VV) : FitFunc(), V(VV) { Np = 2; FitType = "exp"; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V)*exp(-p[1]*t); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/(p[1]*V)*exp(-p[1]*t), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V)*(1.0+p[1]*t)*exp(-p[1]*t) };
  }
  
  virtual ~ExpFunc(){}
};

// Cosh form: y[t] = p[0]^2/(2.0*p[1]*V) * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) )
class CoshFunc : public FitFunc{
  
protected:
  double T;
  double V;
    
public:    
  CoshFunc(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 2; FitType = "cosh"; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V) * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) ); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/(p[1]*V) * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) ), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V) * ( (1.0+p[1]*t)*exp(-p[1]*t) + 
        (1.0+p[1]*(T-t))*exp(-p[1]*(T-t)) ) };
  }
  
  virtual ~CoshFunc(){}
};

#endif