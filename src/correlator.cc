#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "correlator.h"
#include "fitter_controls.h"

std::vector<double> tokenize_line(std::string& s, std::string delim)
{
  std::istringstream iss(s);
  std::string token;
  std::vector<double> tokens;
  
  while(getline(iss,token,' ')){ if(!token.empty()){ tokens.push_back(atof(token.c_str())); } }
  if(tokens.size() != 3){ 
    printf("Error: expected data format is Nx3 matrix with format: (t,real(C(t)),imag(C(t))).\n");
    exit(-1);
  }
  
  return tokens;
}

Correlator::Correlator(fitter_controls& fc, int this_fit)
{
  // Add number of data points in this fit to total
  fc.Ndata += fc.fits[this_fit].t_max - fc.fits[this_fit].t_min + 1;
  
  ////////////////////////////////////////////////////////////////
  // Read correlator data from disk. 
  //
  // This assumes a file containing an Nx2 matrix formatted as
  //    t   C(t)
  // for each trajectory included in the fit.
  // 3-pt functions should be treated as a simultaneous fit of a 
  // series of 2-pt correlators for each source-sink separation.
  ////////////////////////////////////////////////////////////////
  int traj = fc.traj_start;
  fc.Ntraj = ( fc.traj_end - fc.traj_start + fc.traj_inc ) / fc.traj_inc;
  t.resize(fc.Ntraj); 
  C.resize(fc.Ntraj);
  include_in_fit.resize(fc.Ntraj);
  
  for(int i=0; i<fc.Ntraj; ++i)
  {
    // Path to data for this trajectory
    std::stringstream ss;
    ss << fc.fits[this_fit].data_stem << "." << traj;
    
    // Open data file
    std::ifstream in(ss.str());
    if(!in.is_open()){ 
      printf("Error: failed to open file %s\n", ss.str().c_str()); 
      exit(-1); 
    }
    
    // Parse line by line and store result
    std::string line;
    std::vector<double> tokens;
    while(getline(in,line)){ 
      tokens = tokenize_line(line, " ");
      t[i].push_back(tokens[0]);
      C[i].push_back(tokens[1]);
      if((tokens[0] >= fc.fits[this_fit].t_min) && (tokens[0] <= fc.fits[this_fit].t_max)){
        include_in_fit[i].push_back(1);
      } else {
        include_in_fit[i].push_back(0);
      }
    }
    
    // some debug messages
    #if(0)
      printf("\nConstructing correlator for fit %d:\n", this_fit); 
      printf("\n--- Fitting to t = (%d,%d)\n", static_cast<int>(fc.fits[this_fit].t_min), 
        static_cast<int>(fc.fits[this_fit].t_max));
      printf("\n--- Got data for trajectory %d:\n", traj);
      for(unsigned int j=0; j<t[i].size(); ++j){
        printf("%3d\t%1.6e\n", static_cast<int>(t[i][j]), C[i][j]);
      }
    #endif
    
    traj += fc.traj_inc;
  }
  
  // Construct the appropriate FitFunc object for this fit form
  // The available forms are defined in include/fit_functions.h
  double V = static_cast<double>( pow(fc.L,3)*fc.T );
  if(fc.fits[this_fit].fit_type == "linear"){ f = new LinearFunc(); }
  else if(fc.fits[this_fit].fit_type == "exp"){ f = new ExpFunc(V); }
  else if(fc.fits[this_fit].fit_type == "cosh"){ f = new CoshFunc(fc.T,V); }
  else{ 
    printf("Error: unrecognized fit type %s\n", fc.fits[this_fit].fit_type.c_str()); 
    exit(-1); 
  }
  
  fc.Nparams += f->get_Np();
}

Correlator::~Correlator(){ delete f; }