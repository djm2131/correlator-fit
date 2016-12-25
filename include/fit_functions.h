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
  
  virtual double get_mass(std::vector<double>& p) = 0;
    
  virtual double eval(double& t, std::vector<double>& p) = 0;
    
  virtual std::vector<double> eval_derivs(double& t, std::vector<double>& p) = 0;

  virtual ~FitFunc(){}
};

// Constant function: y[t] = p[0]
class ConstFunc : public FitFunc{
    
public:    
  ConstFunc() : FitFunc() { Np = 1; FitType = "constant"; }
  
  double get_mass(std::vector<double>& p){ return p[0]; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return p[0]; 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return {1.0};
  }
  
  virtual ~ConstFunc(){}
};

// Linear function: y[t] = p[0] + p[1]*t
class LinearFunc : public FitFunc{
    
public:    
  LinearFunc() : FitFunc() { Np = 2; FitType = "linear"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
    
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

// Single exponential: y[t] = Z^2/(2*m*V) * e^(-m*t)
//  p[0] = Z
//  p[1] = m
class ExpFunc : public FitFunc{
  
protected:
  double V;
    
public:    
  ExpFunc(double VV) : FitFunc(), V(VV) { Np = 2; FitType = "exp"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V)*exp(-p[1]*t); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V*exp(-p[1]*t), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V)*(1.0+p[1]*t)*exp(-p[1]*t) };
  }
  
  virtual ~ExpFunc(){}
};

// Cosh form for <PP>: y[t] = Z^2/(2*m*V) * ( e^(-m*t) + e^(-m*(T-t)) )
//  p[0] = Z
//  p[1] = m
class CoshFuncPP : public FitFunc{
  
protected:
  double T;
  double V;
    
public:    
  CoshFuncPP(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 2; FitType = "cosh_pp"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V) * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) ); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) ), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V) * ( (1.0+p[1]*t)*exp(-p[1]*t) + 
        (1.0+p[1]*(T-t))*exp(-p[1]*(T-t)) ) };
  }
  
  virtual ~CoshFuncPP(){}
};

// Cosh form for <AP>: y[t] = Z^2/(2*m*V) * ( e^(-m*t) - e^(-m*(T-t)) )
//  p[0] = Z
//  p[1] = m
class CoshFuncAP : public FitFunc{
  
protected:
  double T;
  double V;
    
public:    
  CoshFuncAP(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 2; FitType = "cosh_ap"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V) * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V) * ( (1.0+p[1]*t)*exp(-p[1]*t) - 
        (1.0+p[1]*(T-t))*exp(-p[1]*(T-t)) ) };
  }
  
  virtual ~CoshFuncAP(){}
};

// Cosh form for global fits including decay constants: y[t] = Z*f*sqrt(V)/(2*ZA) * ( e^(-mt) - e^(-m*(T-t)) )
//  p[0] = Z
//  p[1] = m
//  p[2] = ZA
//  p[3] = f
// The first three parameters should be bound to the Z-factor and mass of the corresponding <PP> correlator, and
// the axial current renormalization factor ZA, respectively
class CoshFuncDecay : public FitFunc{
  
protected:
  double T;
  double V;
    
public:    
  CoshFuncDecay(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 4; FitType = "cosh_decay"; }
  
  double get_mass(std::vector<double>& p){ return p[3]; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return 0.5*p[0]*p[3]/p[2] * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { 0.5*p[3]/p[2] * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ), 
      -0.5*p[0]*p[3]/p[2] * ( t*exp(-p[1]*t) - (T-t)*exp(-p[1]*(T-t)) ),
      -0.5*p[0]*p[3]/p[2]/p[2] * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ),
      0.5*p[0]/p[2] * ( exp(-p[1]*t) - exp(-p[1]*(T-t)) ) };
  }
  
  virtual ~CoshFuncDecay(){}
};

// Cosh form for two-pion correlator: y[t] = Z^2/(2*m*V) * ( e^{-m*t} + e^{-m*(T-t)} + C )
//  p[0] = Z
//  p[1] = m
//  p[2] = C
class CoshFuncTwoPion : public FitFunc{

protected:
  double T;
  double V;
  
public:
  CoshFuncTwoPion(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 3; FitType = "cosh_two_pion"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
  
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V) * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) + p[2] ); 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) + p[2] ), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V) * ( (1.0+p[1]*t)*exp(-p[1]*t) + 
        (1.0+p[1]*(T-t))*exp(-p[1]*(T-t)) + p[2] ) ,
        pow(p[0],2.0)/(2.0*p[1]*V) };
  }
  
};

#endif