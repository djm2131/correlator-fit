#ifndef __FITFUNC_H_INCLUDED__
#define __FITFUNC_H_INCLUDED__

#include <string>
#include <cmath>

template <typename T> double sgn(T val){
  return static_cast<double>( (T(0) < val) - (val < T(0)) );
}

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
  
  virtual double thermal_state(double& t, std::vector<double>& p) = 0;
    
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return 0.0; 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return 0.0; 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return 0.0; 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V*exp(-p[1]*t), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V)*(1.0+p[1]*t)*exp(-p[1]*t) };
  }
  
  virtual ~ExpFunc(){}
};

// Two-state exponential: y[t] = Z_1^2/(2*m_1*V) * e^(-m_1*t) + sgn(Z_2)*Z_2^2/(2*m_2*V) * e^(-m_2*t)
// Need to allow for Z_1 and Z_2 to have different relative signs
//  p[0] = Z_1
//  p[1] = m_1
//  p[2] = Z_2
//  p[3] = m_2
class DoubleExpFunc : public FitFunc{
  
protected:
  double V;
    
public:    
  DoubleExpFunc(double VV) : FitFunc(), V(VV) { Np = 4; FitType = "double_exp"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; } // ground state mass
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V)*exp(-p[1]*t) + sgn(p[2])*pow(p[2],2.0)/(2.0*p[3]*V)*exp(-p[3]*t); 
  }
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return 0.0; 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V*exp(-p[1]*t), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V)*(1.0+p[1]*t)*exp(-p[1]*t),
      sgn(p[2])*p[2]/p[3]/V*exp(-p[3]*t), 
      -sgn(p[2])*pow(p[2],2.0)/(2.0*pow(p[3],2.0)*V)*(1.0+p[3]*t)*exp(-p[3]*t) };
  }
  
  virtual ~DoubleExpFunc(){} 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*V) * exp(-p[1]*(T-t)); 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return -pow(p[0],2.0)/(2.0*p[1]*V) * exp(-p[1]*(T-t)); 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    return -0.5*p[0]*p[3]/p[2] * exp(-p[1]*(T-t)); 
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
  
  double thermal_state(double& t, std::vector<double>& p){ 
    check_Np(p);
    // return pow(p[0],2.0)/(2.0*p[1]*V) * ( exp(-p[1]*(T-t)) + p[2] ); 
    return pow(p[0],2.0)/(2.0*p[1]*V)*p[2]; 
  }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/V * ( exp(-p[1]*t) + exp(-p[1]*(T-t)) + p[2] ), 
      -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*V) * ( (1.0+p[1]*t)*exp(-p[1]*t) + 
        (1.0+p[1]*(T-t))*exp(-p[1]*(T-t)) + p[2] ) ,
        pow(p[0],2.0)/(2.0*p[1]*V) };
  }
  
};

// Ratio form for two-pion correlator: R[t+1/2] = Z * ( cosh(dE*(t+1/2-T/2)) + sinh(dE*(t+1/2-T/2))*coth(2*m*(t+1/2-T/2)) )
//  p[0] = Z
//  p[1] = m
//  p[2] = dE
class TwoPionRatioFunc : public FitFunc{

protected:
  double T;
  double V;
  
public:
  TwoPionRatioFunc(double TT, double VV) : FitFunc(), T(TT), V(VV) { Np = 3; FitType = "ratio_two_pion"; }
  
  double get_mass(std::vector<double>& p){ return p[1]; }
  
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    // double tp = t + 0.5*(1.0-T);
    double tp = t - 0.5*T;
    return p[0] * ( cosh(p[2]*tp) + sinh(p[2]*tp)/tanh(2.0*p[1]*tp) ); 
  }
  
  double thermal_state(double& t, std::vector<double>& p){ return 0.0; }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    // double tp = t + 0.5*(1.0-T);
    double tp = t - 0.5*T;
    return { cosh(p[2]*tp) + sinh(p[2]*tp)/tanh(2.0*p[1]*tp), 
      -2.0*p[0]*tp*sinh(p[2]*tp)/pow(sinh(2.0*p[1]*tp),2),
      p[0]*tp*cosh((2.0*p[1]+p[2])*tp)/sinh(2.0*p[1]*tp) };
  }
  
};

// For global fits including ZV measured through EM current: y[t_sep] = 1/ZV * Z^2/(2*m*V) * e^{-m*\Delta t_sep}
//  p[0] = Z
//  p[1] = m
//  p[2] = ZV
class ZVFunc : public FitFunc{
  
protected:
  double T;
  double V;
  double dt; // source-sink separation
    
public:    
  ZVFunc(double TT, double VV, double dtt) : FitFunc(), T(TT), V(VV), dt(dtt) { Np = 3; FitType = "zv"; }
  
  double get_mass(std::vector<double>& p){ return 0.0; }
    
  double eval(double& t, std::vector<double>& p){ 
    check_Np(p);
    return pow(p[0],2.0)/(2.0*p[1]*p[2]*V) * exp(-p[1]*dt); 
  }
  
  double thermal_state(double& t, std::vector<double>& p){ return 0.0; }
  
  std::vector<double> eval_derivs(double& t, std::vector<double>& p){
    check_Np(p);
    return { p[0]/p[1]/p[2]/V * exp(-p[1]*dt), 
             -pow(p[0],2.0)/(2.0*pow(p[1],2.0)*p[2]*V) * (1.0+p[1]*dt)*exp(-p[1]*dt), 
             -pow(p[0],2.0)/(2.0*p[1]*pow(p[2],2.0)*V) * exp(-p[1]*dt) };
  }
  
  virtual ~ZVFunc(){}
};

#endif