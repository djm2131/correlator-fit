#include <cstdio>
#include "xml_parser.h"
#include "fitter_controls.h"

void usage(void){
  printf("Usage: ./fitter <path-to-XML>\n");
  write_xml_template("./template.xml");
}

int main(int argc, char **argv){
  
  if(argc != 2){ usage(); return(0); }

  fitter_controls fit_controls;
  {
    XML_parser XML(argv[1]);
    fit_controls = XML.parse_all();
  }
  
  printf("Got lattice params:\n");
  printf("  L = %d\n", fit_controls.L);
  printf("  T = %d\n", fit_controls.T);
  printf("  traj_start = %d\n", fit_controls.traj_start);
  printf("  traj_end = %d\n", fit_controls.traj_end);
  printf("  traj_inc = %d\n", fit_controls.traj_inc);
  printf("  bin_size = %d\n", fit_controls.bin_size);
  
  printf("\nGot %d fits.\n", static_cast<int>(fit_controls.fits.size()));
  for(unsigned int i=0; i<fit_controls.fits.size(); ++i){
    printf("\nFit %d:\n", i);
    printf("  t_min = %f\n", fit_controls.fits[i].t_min);
    printf("  t_max = %f\n", fit_controls.fits[i].t_max);
    printf("  data_stem = %s\n", fit_controls.fits[i].data_stem.c_str());
    printf("  fit_type = %s\n", fit_controls.fits[i].fit_type.c_str());
  }

  return(0);
}