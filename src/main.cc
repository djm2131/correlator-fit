#include <cstdio>
#include "fitter.h"

void usage(void){
  printf("Usage: ./fitter <path-to-XML>\n");
}

int main(int argc, char **argv){
  
  if(argc != 2){ usage(); return(0); }
  
  Fitter fit(argv[1]);
  fit.do_fit();

  return(0);
}