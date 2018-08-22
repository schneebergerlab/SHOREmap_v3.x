/* this function initialize a SNP struct variable */
#include <string.h>
#include "globals.h"

bool init_SNP(mySNP* iSNP)
{
    (*iSNP).ecotype       = "";
    (*iSNP).chromosome    = "";
    (*iSNP).position      = 0;
    (*iSNP).stype         = "intergenic";
    (*iSNP).gene_id       = "NA";
    (*iSNP).gene_pos      = 0;
    (*iSNP).cds_pos       = 0;
    (*iSNP).codon_pos     = 0;
    (*iSNP).ns_change     = 0; // if !=0, nonsyn effect on aa
    (*iSNP).new_stop      = 0;
    (*iSNP).lost_stop     = 0;
    (*iSNP).splicechange  = 0;
    (*iSNP).ref_base      = "";
    (*iSNP).new_base      = "";
    (*iSNP).ref_aa        = "";
    (*iSNP).new_aa        = "";
    //map<std::string, std::string> domain_change;
    (*iSNP).support       = 0;
    (*iSNP).concordance   = 0;
    (*iSNP).quality       = 0;
    (*iSNP).peak_distance = 999999999;
    (*iSNP).marker_ratio  = 0;
    (*iSNP).source        = "IlluminaGA2";
    
    return true;
}
