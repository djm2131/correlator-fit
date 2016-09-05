#ifndef __XML_PARSER_H_INCLUDED__
#define __XML_PARSER_H_INCLUDED__

#include <cstdio>
#include <string>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include "fitter_controls.h"

void write_xml_template(std::string fout){
  printf("Writing XML template to %s...\n", fout.c_str());
  FILE* fp = fopen(fout.c_str(), "w");
  fprintf(fp, "<?xml version='1.0' encoding='us-ascii'?>\n");
  fprintf(fp, "<fit>\n");
  fprintf(fp, "\t<Lattice>\n");
  fprintf(fp, "\t\t<L>16</L>\n");
  fprintf(fp, "\t\t<T>32</T>\n");
  fprintf(fp, "\t</Lattice>\n");
  fprintf(fp, "\t<correlator>\n");
  fprintf(fp, "\t\t<t_min>2</t_min>\n");
  fprintf(fp, "\t\t<t_max>8</t_max>\n");
  fprintf(fp, "\t\t<data_path>PATH-TO-DATA</data_path>\n");
  fprintf(fp, "\t\t<fit_func>linear</fit_func>\n");
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
  
public:
  XML_parser(){}
  
  XML_parser(std::string path) : path_to_doc(path){
    doc = xmlParseFile(path_to_doc.c_str());
    if(doc == NULL){ parse_error(path_to_doc); }
  }
  
  fitter_controls parse_all(void){
    
    fitter_controls fc;
    xmlXPathContext* xpathCtx;
    xmlXPathObject* xpathObj;
    
    xpathCtx = xmlXPathNewContext(doc);
    if(xpathCtx == NULL){ printf("Error: failed to create new XPath context\n"); exit(-1); }
    
    xpathObj = xmlXPathEvalExpression((xmlChar*) "/fit/Lattice/L", xpathCtx);
    fc.L = atof((char*) xmlNodeListGetString(doc, xpathObj->nodesetval->nodeTab[0]->xmlChildrenNode, 1));
    
    xpathObj = xmlXPathEvalExpression((xmlChar*) "/fit/Lattice/T", xpathCtx);
    fc.T = atof((char*) xmlNodeListGetString(doc, xpathObj->nodesetval->nodeTab[0]->xmlChildrenNode, 1));
    
    xmlXPathFreeObject(xpathObj);
    xmlXPathFreeContext(xpathCtx);
    
    return fc;
  }
  
  ~XML_parser(){ xmlFreeDoc(doc); }
};

#endif