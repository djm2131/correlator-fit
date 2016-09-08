#ifndef __FITFUNC_H_INCLUDED__
#define __FITFUNC_H_INCLUDED__

#include <string>

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
    // std::vector<double> d = {1.0, t};
    // return d;
    return {1.0, t};
  }
  
  virtual ~LinearFunc(){}
};

#endif