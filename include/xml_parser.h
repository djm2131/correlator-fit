#ifndef __XML_PARSER_H_INCLUDED__
#define __XML_PARSER_H_INCLUDED__

#include <cstdio>
#include <string>
#include <libxml/xpath.h>
#include "fitter_controls.h"

class XML_parser{
  
private:
  std::string path_to_doc;
  xmlDoc* doc;
  xmlXPathContext* xpc;
  xmlXPathObject* xpo;
  
public:

  void write_xml_template(std::string fout)
  {
    printf("Writing XML template to %s...\n", fout.c_str());
    FILE* fp = fopen(fout.c_str(), "w");
    fprintf(fp, "<?xml version='1.0' encoding='us-ascii'?>\n");
    fprintf(fp, "<fit>\n");
    fprintf(fp, "\t<Lattice>\n");
    fprintf(fp, "\t\t<L>16</L>\n");
    fprintf(fp, "\t\t<T>32</T>\n");
    fprintf(fp, "\t\t<traj_start>0</traj_start>\n");
    fprintf(fp, "\t\t<traj_end>0</traj_end>\n");
    fprintf(fp, "\t\t<traj_inc>1</traj_inc>\n");
    fprintf(fp, "\t\t<bin_size>1</bin_size>\n");
    fprintf(fp, "\t</Lattice>\n");
    fprintf(fp, "\t<correlator>\n");
    fprintf(fp, "\t\t<t_min>2</t_min>\n");
    fprintf(fp, "\t\t<t_max>8</t_max>\n");
    fprintf(fp, "\t\t<data_stem>./data0</data_stem>\n");
    fprintf(fp, "\t\t<fit_type>linear</fit_type>\n");
    fprintf(fp, "\t</correlator>\n");
    fprintf(fp, "\t<correlator>\n");
    fprintf(fp, "\t\t<t_min>3</t_min>\n");
    fprintf(fp, "\t\t<t_max>7</t_max>\n");
    fprintf(fp, "\t\t<data_stem>./data1</data_stem>\n");
    fprintf(fp, "\t\t<fit_type>linear</fit_type>\n");
    fprintf(fp, "\t</correlator>\n");
    fprintf(fp, "</fit>");
    fclose(fp);
  }

  void parse_error(std::string path)
  {
    printf("Error: could not parse file %s\n", path.c_str());
    write_xml_template("./template.xml");
    exit(-1);
  }

  XML_parser(std::string path) : path_to_doc(path)
  {
    doc = xmlParseFile(path_to_doc.c_str());
    if(doc == NULL){ parse_error(path_to_doc); }
    xpc = xmlXPathNewContext(doc);
    if(xpc == NULL){ printf("Error: failed to create new XPath context\n"); exit(-1); }
  }
  
  int get_Nnodes(const char* path)
  {
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return xmlChildElementCount(xpo->nodesetval->nodeTab[0]);
  }
  
  double parse_numeric(const char* path)
  {
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return atof( reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->xmlChildrenNode, 1) ) );
  }
  
  std::string parse_text(const char* path)
  {
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->xmlChildrenNode, 1) );
  }
  
  void print_params(fitter_controls& fit_controls)
  {
    printf("\nLattice parameters:\n");
    printf("  L = %d\n", fit_controls.L);
    printf("  T = %d\n", fit_controls.T);
    printf("  traj_start = %d\n", fit_controls.traj_start);
    printf("  traj_end = %d\n", fit_controls.traj_end);
    printf("  traj_inc = %d\n", fit_controls.traj_inc);
    printf("  Ntraj = %d\n", fit_controls.Ntraj);
    printf("  bin_size = %d\n", fit_controls.bin_size);
  
    printf("\nFound %d fit(s).\n", static_cast<int>(fit_controls.fits.size()));
    for(unsigned int i=0; i<fit_controls.fits.size(); ++i){
      printf("\nFit %d:\n", i);
      printf("  t_min = %d\n", static_cast<int>(fit_controls.fits[i].t_min));
      printf("  t_max = %d\n", static_cast<int>(fit_controls.fits[i].t_max));
      printf("  data_stem = %s\n", fit_controls.fits[i].data_stem.c_str());
      printf("  fit_type = %s\n", fit_controls.fits[i].fit_type.c_str());
      printf("  initial guesses:\n");
      for(unsigned int j=0; j<fit_controls.p0[i].size(); ++j){
        printf("    p0[%d] = %1.4e\n", j, fit_controls.p0[i][j]);
      }
    }
  }
  
  fitter_controls parse_all(void)
  {
    fitter_controls fit_controls;
    
    // Parse the lattice parameters
    fit_controls.L = parse_numeric("/fit/Lattice/L");
    fit_controls.T = parse_numeric("/fit/Lattice/T");
    fit_controls.traj_start = parse_numeric("/fit/Lattice/traj_start");
    fit_controls.traj_end = parse_numeric("/fit/Lattice/traj_end");
    fit_controls.traj_inc = parse_numeric("/fit/Lattice/traj_inc");
    fit_controls.bin_size = parse_numeric("/fit/Lattice/bin_size");
    fit_controls.Ntraj = ( fit_controls.traj_end - fit_controls.traj_start ) / fit_controls.traj_inc + 1;
    
    // Fill fit_controls.fits (fit_controls.p0) with the details (initial guesses) of each fit
    int Nfits = get_Nnodes("/fit") - 1; // don't count /fit/Lattice
    fit_controls.p0.resize(Nfits);
    fit_controls.fits.resize(Nfits);
    for(int i=0; i<Nfits; ++i){
      
      // fit parameters
      char path[100];
      sprintf(path, "/fit/correlator[%d]/t_min", i+1);     fit_controls.fits[i].t_min = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/t_max", i+1);     fit_controls.fits[i].t_max = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/data_stem", i+1); fit_controls.fits[i].data_stem = parse_text(path);
      sprintf(path, "/fit/correlator[%d]/fit_type", i+1);  fit_controls.fits[i].fit_type = parse_text(path);
      
      // initial guesses
      sprintf(path, "/fit/correlator[%d]/parameter_guesses", i+1);
      int Np = get_Nnodes(path);
      fit_controls.p0[i].resize(Np);
      for(int j=0; j<Np; ++j){ 
        sprintf(path, "/fit/correlator[%d]/parameter_guesses/p%d", i+1, j); 
        fit_controls.p0[i][j] = parse_numeric(path); 
      }
    }
    
    print_params(fit_controls);
    
    return fit_controls;
  }
  
  ~XML_parser()
  { 
    xmlXPathFreeObject(xpo);
    xmlXPathFreeContext(xpc);
    xmlFreeDoc(doc);
    xmlCleanupParser();
  }
};

#endif