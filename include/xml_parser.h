#ifndef __XML_PARSER_H_INCLUDED__
#define __XML_PARSER_H_INCLUDED__

#include <cstdio>
#include <string>
#include <libxml/xpath.h>
#include "fitter_controls.h"

void write_xml_template(std::string fout){
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

void parse_error(std::string path){
  printf("Error: could not parse file %s\n", path.c_str());
  write_xml_template("./template.xml");
  exit(-1);
}

class XML_parser{
  
private:
  std::string path_to_doc;
  xmlDoc* doc;
  xmlXPathContext* xpc;
  xmlXPathObject* xpo;
  
public:
  XML_parser(){}
  
  XML_parser(std::string path) : path_to_doc(path){
    doc = xmlParseFile(path_to_doc.c_str());
    if(doc == NULL){ parse_error(path_to_doc); }
    xpc = xmlXPathNewContext(doc);
    if(xpc == NULL){ printf("Error: failed to create new XPath context\n"); exit(-1); }
  }
  
  int get_Nfits(void){
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>("/fit") );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return xmlChildElementCount(xpo->nodesetval->nodeTab[0]) - 1;
  }
  
  double parse_numeric(const char* path){
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return atof( reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->xmlChildrenNode, 1) ) );
  }
  
  std::string parse_text(const char* path){
    xmlChar* xpath = const_cast<xmlChar*>( reinterpret_cast<const xmlChar*>(path) );
    xpo = xmlXPathEvalExpression(xpath, xpc);
    return reinterpret_cast<char*>( xmlNodeListGetString(doc, xpo->nodesetval->nodeTab[0]->xmlChildrenNode, 1) );
  }
  
  fitter_controls parse_all(void){
    
    fitter_controls fc;
    
    // Parse the lattice parameters
    fc.L = parse_numeric("/fit/Lattice/L");
    fc.T = parse_numeric("/fit/Lattice/T");
    fc.traj_start = parse_numeric("/fit/Lattice/traj_start");
    fc.traj_end = parse_numeric("/fit/Lattice/traj_end");
    fc.traj_inc = parse_numeric("/fit/Lattice/traj_inc");
    fc.bin_size = parse_numeric("/fit/Lattice/bin_size");
    
    // Fill fc.fits with the details of each fit
    fc.fits.resize(get_Nfits());
    for(unsigned int i=0; i<fc.fits.size(); ++i){
      char path[100];
      sprintf(path, "/fit/correlator[%d]/t_min", i+1);     fc.fits[i].t_min = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/t_max", i+1);     fc.fits[i].t_max = parse_numeric(path);
      sprintf(path, "/fit/correlator[%d]/data_stem", i+1); fc.fits[i].data_stem = parse_text(path);
      sprintf(path, "/fit/correlator[%d]/fit_type", i+1);  fc.fits[i].fit_type = parse_text(path);
    }
    
    return fc;
  }
  
  ~XML_parser(){ 
    xmlXPathFreeObject(xpo);
    xmlXPathFreeContext(xpc);
    xmlFreeDoc(doc);
  }
};

#endif