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
  
  printf("Got %d^3x%d lattice!\n", fit_controls.L, fit_controls.T);

  return(0);
}