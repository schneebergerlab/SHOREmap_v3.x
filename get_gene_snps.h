/* this function set up the CDS sequence for a specific gene for checking protein changes */

#include "globals.h"

/* singleGENEann = gene_ann[gene_name]; snp_list from get_SNPlist(...) */
void get_gene_snps(std::string chromosome, 
                   std::string gene_name,
                   QUARTET singleGENEann, 
                   map<unsigned long, QUARTET> gene_coding_ann,
                   std::string isoform,
                   multimap<unsigned long, mySNP> snp_list,
                   GeneSNPlist* gene,
                   FILE* fplog);
