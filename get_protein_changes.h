/* this function calculates CDS&proteins changes of all SNPs in SNP_lists and fills in SNP object */
/* date: 2013-04-15                                                                               */

#include "globals.h"
void get_protein_changes(GeneSNPlist* gene, multimap<unsigned long, mySNP>* tSNPlist);
