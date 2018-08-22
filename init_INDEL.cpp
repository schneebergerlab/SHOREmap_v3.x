#include <string.h>
#include "globals.h"

bool init_INDEL(myINDEL* iINDEL)
{
    (*iINDEL).ecotype       = "";
    (*iINDEL).chromosome    = "";
    (*iINDEL).begin         =  0;
    (*iINDEL).end           =  0;
    (*iINDEL).seq           = "";
    (*iINDEL).stype         = "intergenic";
    
    (*iINDEL).gene_id       = "NA";
    (*iINDEL).gene_pos      = 0;
    (*iINDEL).cds_pos       = 0;
    (*iINDEL).codon_pos     = 0;
    (*iINDEL).new_stop      = 0;
    (*iINDEL).lost_stop     = 0;
    (*iINDEL).splicechange  = 0;
    //map<std::string, std::string> domain_change;
    
    (*iINDEL).support       = 0;
    (*iINDEL).concordance   = 0;
    (*iINDEL).quality       = 0;
    (*iINDEL).peak_distance = 999999999;
    (*iINDEL).marker_ratio  = 0;
    (*iINDEL).source        = "IlluminaGA2";
    
    return true;
}

/*

    std::string   ecotype;
    std::string   chromosome;
    unsigned long begin;
    unsigned long end;
    stf::string   seq;
    std::string   stype;
    std::string   gene_id;
    unsigned long gene_pos;  // 
    unsigned long cds_pos;   // 
    unsigned long codon_pos; //
    unsigned long ns_change;
    unsigned long new_stop;
    unsigned long lost_stop;
    unsigned long splicechange;
    //map<std::string, std::string> domain_change;
    unsigned long support;
    unsigned long support;
    double        concordance;
    double        quality;
    unsigned long peak_distance;
    double        marker_ratio;
    std::string   source;

*/
