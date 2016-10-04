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
    fprintf(fp, "\t<LM_fitter>\n");
    fprintf(fp, "\t\t<numerical_derivs>0</numerical_derivs>\n");
    fprintf(fp, "\t\t<max_iter>1000</max_iter>\n");
    fprintf(fp, "\t\t<xtol>1.0e-08</xtol>\n");
    fprintf(fp, "\t\t<gtol>1.0e-08</gtol>\n");
    fprintf(fp, "\t\t<ftol>0.0</ftol>\n");
    fprintf(fp, "\t</LM_fitter>\n");
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
    fprintf(fp, "\t\t<parameter_guesses>\n");
    fprintf(fp, "\t\t\t<p0 name=\"a\">1.0</p0>\n");
    fprintf(fp, "\t\t\t<p1 name=\"b\">-1.0</p1>\n");
    fprintf(fp, "\t\t</parameter_guesses>\n");
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
  
  std::string parse_attribute(const char* path)
  {
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->properties->children, 1) );
  }
  
  void print_params(fitter_controls& fit_controls)
  {
    printf("\nLevenberg-Marquardt fitter parameters:\n");
    printf("  numerical_derivs = %d\n", fit_controls.numerical_derivs);
    printf("  max_iter = %d\n", fit_controls.max_iter);
    printf("  xtol = %e\n", fit_controls.xtol);
    printf("  gtol = %e\n", fit_controls.gtol);
    printf("  ftol = %e\n", fit_controls.ftol);
    
    printf("\nLattice parameters:\n");
    printf("  L = %d\n", fit_controls.L);
    printf("  T = %d\n", fit_controls.T);
    printf("  traj_start = %d\n", fit_controls.traj_start);
    printf("  traj_end = %d\n", fit_controls.traj_end);
    printf("  traj_inc = %d\n", fit_controls.traj_inc);
    printf("  Ntraj = %d\n", fit_controls.Ntraj);
    printf("  bin_size = %d\n", fit_controls.bin_size);
  
    printf("\nFound %d correlator(s).\n", static_cast<int>(fit_controls.fits.size()));
    int p_idx(0);
    for(unsigned int i=0; i<fit_controls.fits.size(); ++i){
      printf("\nCorrelator %d:\n", i);
      printf("  t_min = %d\n", static_cast<int>(fit_controls.fits[i].t_min));
      printf("  t_max = %d\n", static_cast<int>(fit_controls.fits[i].t_max));
      printf("  data_stem = %s\n", fit_controls.fits[i].data_stem.c_str());
      printf("  fit_type = %s\n", fit_controls.fits[i].fit_type.c_str());
      printf("  initial guesses:\n");
      for(unsigned int j=0; j<fit_controls.p0[i].size(); ++j){
        printf("    parameter %d: %s = %1.4e\n", j, fit_controls.p_names[p_idx].c_str(), fit_controls.p0[i][j]);
        ++p_idx;
      }
    }
  }
  
  fitter_controls parse_all(void)
  {
    fitter_controls fit_controls;
    
    // Parse the Levenberg-Marquardt parameters
    fit_controls.numerical_derivs = parse_numeric("/fit/LM_fitter/numerical_derivs");
    fit_controls.max_iter = parse_numeric("/fit/LM_fitter/max_iter");
    fit_controls.xtol = parse_numeric("/fit/LM_fitter/xtol");
    fit_controls.gtol = parse_numeric("/fit/LM_fitter/gtol");
    fit_controls.ftol = parse_numeric("/fit/LM_fitter/ftol");
    
    // Parse the lattice parameters
    fit_controls.L = parse_numeric("/fit/Lattice/L");
    fit_controls.T = parse_numeric("/fit/Lattice/T");
    fit_controls.traj_start = parse_numeric("/fit/Lattice/traj_start");
    fit_controls.traj_end = parse_numeric("/fit/Lattice/traj_end");
    fit_controls.traj_inc = parse_numeric("/fit/Lattice/traj_inc");
    fit_controls.bin_size = parse_numeric("/fit/Lattice/bin_size");
    fit_controls.Ntraj = ( fit_controls.traj_end - fit_controls.traj_start ) / fit_controls.traj_inc + 1;
    
    // Fill fit_controls.fits (fit_controls.p0) with the details (initial guesses) of each fit
    int Nfits = get_Nnodes("/fit") - 3; // don't count LM, lattice params, or constraints
    fit_controls.p0.resize(Nfits);
    fit_controls.fits.resize(Nfits);
    char path[100];
    for(int i=0; i<Nfits; ++i)
    {
      // fit parameters
      sprintf(path, "/fit/correlator[%d]/t_min", i+1);     fit_controls.fits[i].t_min = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/t_max", i+1);     fit_controls.fits[i].t_max = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/data_stem", i+1); fit_controls.fits[i].data_stem = parse_text(path);
      sprintf(path, "/fit/correlator[%d]/fit_type", i+1);  fit_controls.fits[i].fit_type = parse_text(path);
      
      // initial guesses and parameter names
      sprintf(path, "/fit/correlator[%d]/parameter_guesses", i+1);
      int Np = get_Nnodes(path);
      fit_controls.p0[i].resize(Np);
      for(int j=0; j<Np; ++j){ 
        sprintf(path, "/fit/correlator[%d]/parameter_guesses/p%d", i+1, j); 
        fit_controls.p0[i][j] = parse_numeric(path); 
        fit_controls.p_names.push_back(parse_attribute(path));
      }
    }
    fit_controls.p_names.push_back("chisq/dof");
    
    print_params(fit_controls);
    
    // Parse parameter bindings (if any are specified)
    sprintf(path, "/fit/Constraints");
    int Nc = get_Nnodes(path);
    if(Nc > 0){
      fit_controls.constrained_fit = true;
      for(int i=0; i<Nc; ++i){ 
        std::vector<int> this_constraint(4);
        sprintf(path, "/fit/Constraints/param_bind[%d]/corr_1", i+1);  this_constraint[0] = parse_numeric(path);
        sprintf(path, "/fit/Constraints/param_bind[%d]/param_1", i+1); this_constraint[1] = parse_numeric(path);
        sprintf(path, "/fit/Constraints/param_bind[%d]/corr_2", i+1);  this_constraint[2] = parse_numeric(path);
        sprintf(path, "/fit/Constraints/param_bind[%d]/param_2", i+1); this_constraint[3] = parse_numeric(path);
        fit_controls.p_bindings.push_back(this_constraint); 
        printf("\nBinding correlator %d parameter %d to correlator %d parameter %d.", 
          this_constraint[0], this_constraint[1], this_constraint[2], this_constraint[3]);
      }
      printf("\n");
    } else {
      fit_controls.constrained_fit = false;
      printf("\nPerforming unconstrained fit.\n");
    }
    
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