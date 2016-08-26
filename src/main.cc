#include <cstdio>
#include <string>
#include <libxml/parser.h>
#include <libxml/tree.h>

void write_xml_template(std::string fout)
{
    printf("Writing XML template to %s...\n", fout.c_str());
    
    FILE* fp = fopen(fout.c_str(), "w");
    fprintf(fp, "<fit>\n");
    fprintf(fp, "\t<Lattice>\n");
    fprintf(fp, "\t\t<L>16<\\L>\n");
    fprintf(fp, "\t\t<T>32<\\T>\n");
    fprintf(fp, "\t<\\Lattice>\n");
    fprintf(fp, "\t<correlator>\n");
    fprintf(fp, "\t\t<t_min>2<\\t_min>\n");
    fprintf(fp, "\t\t<t_max>8<\\t_max>\n");
    fprintf(fp, "\t\t<data_path>PATH-TO-DATA</data_path>\n");
    fprintf(fp, "\t\t<fit_func>linear<\\fit_func>\n");
    fprintf(fp, "\t<\\correlator>\n");
    fprintf(fp, "<\\fit>");
    fclose(fp);
}

int main(int argc, char **argv)
{
    if(argc != 2){  
        printf("Usage: ./fitter <path-to-XML>\n");
        write_xml_template("./template.xml");
        return(0);
    }

    xmlDoc *doc = NULL;
    doc = xmlReadFile(argv[1], NULL, 0);
    if(doc == NULL){ printf("Error: could not parse file %s\n", argv[1]); }

    xmlNode *root_element = xmlDocGetRootElement(doc);

    xmlFreeDoc(doc);

    xmlCleanupParser();

    return(0);
}