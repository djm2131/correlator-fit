/**********************************************************************
    
    CorrelatorFit (www.github.com/djm2131/correlator-fit)
    
    Source file: include/xml_parser.h
    
    Author: David Murphy (djmurphy@mit.edu)
    
    Copyright (c) 2019 MIT. All rights reserved.

**********************************************************************/

#ifndef __XML_PARSER_H_INCLUDED__
#define __XML_PARSER_H_INCLUDED__

#include <cstdlib>
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

    void write_xml_template(const std::string& fout) const
    {
      printf("Writing XML template to %s...\n", fout.c_str());
      FILE* fp = fopen(fout.c_str(), "w");
      fprintf(fp, "<?xml version='1.0' encoding='us-ascii'?>\n");
      fprintf(fp, "<fit>\n");
      fprintf(fp, "\t<fitter>\n");
      fprintf(fp, "\t\t<algorithm>NM</algorithm>\n");
      fprintf(fp, "\t\t<numerical_derivs>0</numerical_derivs>\n");
      fprintf(fp, "\t\t<max_iter>5000</max_iter>\n");
      fprintf(fp, "\t\t<xtol>1.0e-08</xtol>\n");
      fprintf(fp, "\t\t<gtol>1.0e-08</gtol>\n");
      fprintf(fp, "\t\t<ftol>0.0</ftol>\n");
      fprintf(fp, "\t\t<correlated_fits>0</correlated_fits>\n");
      fprintf(fp, "\t\t<svd_cut>1.0e-15</svd_cut>\n");
      fprintf(fp, "\t</fitter>\n");
      fprintf(fp, "\t<Lattice>\n");
      fprintf(fp, "\t\t<L>16</L>\n");
      fprintf(fp, "\t\t<T>32</T>\n");
      fprintf(fp, "\t\t<traj_start>0</traj_start>\n");
      fprintf(fp, "\t\t<traj_end>0</traj_end>\n");
      fprintf(fp, "\t\t<traj_inc>1</traj_inc>\n");
      fprintf(fp, "\t\t<bin_size>1</bin_size>\n");
      fprintf(fp, "\t</Lattice>\n");
      fprintf(fp, "\t<correlator>\n");
      fprintf(fp, "\t\t<resample>1</resample>\n");
      fprintf(fp, "\t\t<t_min>2</t_min>\n");
      fprintf(fp, "\t\t<t_max>8</t_max>\n");
      fprintf(fp, "\t\t<t_sep>0</t_sep>\n");
      fprintf(fp, "\t\t<data_stem>data</data_stem>\n");
      fprintf(fp, "\t\t<cov_matrix_stem>mcov.dat</cov_matrix_stem>");
      fprintf(fp, "\t\t<fit_type>linear</fit_type>\n");
      fprintf(fp, "\t\t<parameter_guesses>\n");
      fprintf(fp, "\t\t\t<p0 name=\"a\">1.0</p0>\n");
      fprintf(fp, "\t\t\t<p1 name=\"b\">-1.0</p1>\n");
      fprintf(fp, "\t\t</parameter_guesses>\n");
      fprintf(fp, "\t\t<eff_mass>\n");
      fprintf(fp, "\t\t\t<compute_eff_mass>1</compute_eff_mass>");
      fprintf(fp, "\t\t\t<subtract_thermal_state>0</subtract_thermal_state>");
      fprintf(fp, "\t\t\t<eff_mass_type>log</eff_mass_type>");
      fprintf(fp, "\t\t\t<eff_mass_stem>eff_mass.dat</eff_mass_stem>");
      fprintf(fp, "\t\t</eff_mass>\n");
      fprintf(fp, "\t</correlator>\n");
      fprintf(fp, "\t<Constraints>\n");
      fprintf(fp, "\t</Constraints>\n");
      fprintf(fp, "\t<save_jacks>1</save_jacks>\n");
      fprintf(fp, "\t<jacks_dir>.</jacks_dir>\n");
      fprintf(fp, "\t<save_chi2pdof>1</save_chi2pdof>\n");
      fprintf(fp, "\t<chi2pdof_stem>chi2pdof</chi2pdof_stem>\n");
      fprintf(fp, "</fit>");
      fclose(fp);
    }

    void parse_error(const std::string& path) const
    {
      printf("Error: could not parse file %s\n", path.c_str());
      write_xml_template("./template.xml");
      exit(-1);
    }

    explicit XML_parser(const std::string& path) : path_to_doc(path), doc(nullptr), xpc(nullptr), xpo(nullptr)
    {
      doc = xmlParseFile(path_to_doc.c_str());
      if(doc == nullptr){ parse_error(path_to_doc); }
      xpc = xmlXPathNewContext(doc);
      if(xpc == nullptr){ printf("Error: failed to create new XPath context\n"); exit(-1); }
    }

    size_t get_Nnodes(const char* path)
    {
      xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
      xpo = xmlXPathEvalExpression(xpath, xpc);
      return xmlChildElementCount(xpo->nodesetval->nodeTab[0]);
    }

    double parse_numeric(const char* path)
    {
      xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
      xpo = xmlXPathEvalExpression(xpath, xpc);
      return strtod( reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->xmlChildrenNode, 1) ), nullptr );
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

    void print_params(const fitter_controls& fit_controls) const
    {
      printf("\nLevenberg-Marquardt fitter parameters:\n");
      printf("  numerical_derivs = %d\n", fit_controls.numerical_derivs);
      printf("  max_iter = %d\n", fit_controls.max_iter);
      printf("  xtol = %e\n", fit_controls.xtol);
      printf("  gtol = %e\n", fit_controls.gtol);
      printf("  ftol = %e\n", fit_controls.ftol);
      printf("  correlated_fits = %d\n", static_cast<int>(fit_controls.correlated_fits));
      printf("  svd_cut = %g\n", fit_controls.svd_cut);

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
        printf("  cov_matrix_stem = %s\n", fit_controls.fits[i].cov_matrix_stem.c_str());
        if(fit_controls.fits[i].do_eff_mass){
          printf("  eff_mass_stem = %s\n", fit_controls.fits[i].eff_mass_stem.c_str());
          printf("  eff_mass_type = %s\n", fit_controls.fits[i].eff_mass_type.c_str());
        }
        printf("  initial guesses:\n");
        for(unsigned int j=0; j<fit_controls.p0[i].size(); ++j){
          printf("    parameter %d: %s = %1.4e\n", j, fit_controls.p_names[p_idx].c_str(), fit_controls.p0[i][j]);
          ++p_idx;
        }
      }
    }

    fitter_controls parse_all()
    {
      fitter_controls fit_controls;

      // Parse the Levenberg-Marquardt parameters
      fit_controls.algorithm        = parse_text("/fit/fitter/algorithm");
      fit_controls.numerical_derivs = static_cast<int>( parse_numeric("/fit/fitter/numerical_derivs") );
      fit_controls.max_iter         = static_cast<int>( parse_numeric("/fit/fitter/max_iter") );
      fit_controls.xtol             = parse_numeric("/fit/fitter/xtol");
      fit_controls.gtol             = parse_numeric("/fit/fitter/gtol");
      fit_controls.ftol             = parse_numeric("/fit/fitter/ftol");
      fit_controls.correlated_fits  = static_cast<bool>( parse_numeric("/fit/fitter/correlated_fits") );
      fit_controls.svd_cut          = parse_numeric("/fit/fitter/svd_cut");

      if((fit_controls.algorithm == "LM") && (fit_controls.correlated_fits == true)){
        printf("\nError: Levenberg-Marquardt fitter does not currently support correlated fits.\n");
        exit(-1);
      }

      // Parse the lattice parameters
      fit_controls.L          = static_cast<int>( parse_numeric("/fit/Lattice/L") );
      fit_controls.T          = static_cast<int>( parse_numeric("/fit/Lattice/T") );
      fit_controls.traj_start = static_cast<int>( parse_numeric("/fit/Lattice/traj_start") );
      fit_controls.traj_end   = static_cast<int>( parse_numeric("/fit/Lattice/traj_end") );
      fit_controls.traj_inc   = static_cast<int>( parse_numeric("/fit/Lattice/traj_inc") );
      fit_controls.Ntraj      = ( fit_controls.traj_end - fit_controls.traj_start ) / fit_controls.traj_inc + 1;
      fit_controls.bin_size   = 1;

      // Fill fit_controls.fits (fit_controls.p0) with the details (initial guesses) of each fit
      size_t Nfits = get_Nnodes("/fit") - 7; // TODO: this is not a great design....
      fit_controls.p0.resize(Nfits);
      fit_controls.fits.resize(Nfits);
      char path[100];
      for(int i=0; i<Nfits; ++i)
      {
        // fit parameters
        sprintf(path, "/fit/correlator[%d]/resample", i+1);        fit_controls.fits[i].resample        = static_cast<bool>( parse_numeric(path) );
        sprintf(path, "/fit/correlator[%d]/t_min", i+1);           fit_controls.fits[i].t_min           = parse_numeric(path);
        sprintf(path, "/fit/correlator[%d]/t_max", i+1);           fit_controls.fits[i].t_max           = parse_numeric(path);
        sprintf(path, "/fit/correlator[%d]/t_sep", i+1);           fit_controls.fits[i].t_sep           = parse_numeric(path);
        sprintf(path, "/fit/correlator[%d]/data_stem", i+1);       fit_controls.fits[i].data_stem       = parse_text(path);
        sprintf(path, "/fit/correlator[%d]/cov_matrix_stem", i+1); fit_controls.fits[i].cov_matrix_stem = parse_text(path);
        sprintf(path, "/fit/correlator[%d]/fit_type", i+1);        fit_controls.fits[i].fit_type        = parse_text(path);

        // effective mass controls
        sprintf(path, "/fit/correlator[%d]/eff_mass/compute_eff_mass", i+1);       fit_controls.fits[i].do_eff_mass   = static_cast<bool>( parse_numeric(path) );
        sprintf(path, "/fit/correlator[%d]/eff_mass/subtract_thermal_state", i+1); fit_controls.fits[i].subtract_ts   = static_cast<bool>( parse_numeric(path) );
        sprintf(path, "/fit/correlator[%d]/eff_mass/eff_mass_type", i+1);          fit_controls.fits[i].eff_mass_type = parse_text(path);
        sprintf(path, "/fit/correlator[%d]/eff_mass/eff_mass_stem", i+1);          fit_controls.fits[i].eff_mass_stem = parse_text(path);

        // initial guesses and parameter names
        sprintf(path, "/fit/correlator[%d]/parameter_guesses", i+1);
        size_t Np = get_Nnodes(path);
        fit_controls.p0[i].resize(Np);
        for(int j=0; j<Np; ++j){
          sprintf(path, "/fit/correlator[%d]/parameter_guesses/p%d", i+1, j);
          fit_controls.p0[i][j] = parse_numeric(path);
          fit_controls.p_names.push_back(parse_attribute(path));
        }
      }
      fit_controls.p_names.emplace_back("chisq/dof");

      print_params(fit_controls);

      // Parse parameter bindings (if any are specified)
      sprintf(path, "/fit/Constraints");
      size_t Nc = get_Nnodes(path);
      if(Nc > 0) {
        fit_controls.constrained_fit = true;
        for(int i=0; i<Nc; ++i){
          std::vector<int> this_constraint(4);
          sprintf(path, "/fit/Constraints/param_bind[%d]/corr_1", i+1);  this_constraint[0] = static_cast<int>( parse_numeric(path) );
          sprintf(path, "/fit/Constraints/param_bind[%d]/param_1", i+1); this_constraint[1] = static_cast<int>( parse_numeric(path) );
          sprintf(path, "/fit/Constraints/param_bind[%d]/corr_2", i+1);  this_constraint[2] = static_cast<int>( parse_numeric(path) );
          sprintf(path, "/fit/Constraints/param_bind[%d]/param_2", i+1); this_constraint[3] = static_cast<int>( parse_numeric(path) );
          fit_controls.p_bindings.push_back(this_constraint);
          printf("\nBinding correlator %d parameter %d to correlator %d parameter %d.",
                 this_constraint[0], this_constraint[1], this_constraint[2], this_constraint[3]);
        }
        printf("\n");
      } else {
        fit_controls.constrained_fit = false;
        printf("\nPerforming unconstrained fit.\n");
      }

      // Parse output options
      fit_controls.save_jacks    = static_cast<bool>( parse_numeric("/fit/save_jacks") );
      fit_controls.jacks_dir     = parse_text("/fit/jacks_dir");
      fit_controls.save_chi2pdof = static_cast<bool>( parse_numeric("/fit/save_chi2pdof") );
      fit_controls.chi2pdof_stem = parse_text("/fit/chi2pdof_stem");

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
