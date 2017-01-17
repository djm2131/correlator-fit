#include <cstdio>
#include "fitter.h"
#include "fit_results.h"

void usage(void){ printf("Usage: ./fitter <path-to-XML>\n"); }

int main(int argc, char **argv)
{  
  if(argc != 2){ usage(); return(0); }
  
  // Construct fitter and do fit
  Fitter fit(argv[1]);
  fit_results fr = fit.do_fit();
  fit.print_results(fr);
  
  // Save fit parameters (for global fits)
  fit.save_jacks(fr);
  
  // Compute and save effective masses
  fit.compute_eff_mass(fr);
  fit.save_eff_mass(fr);

  return(0);
}