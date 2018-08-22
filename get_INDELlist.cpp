/* this function reads info of INDELs from given insertion or deletion file */

#include <map>
#include <stdlib.h>
#include <stdio.h>
#include "init_INDEL.h"

void get_INDELlist(char* indelfile, multimap<unsigned long, myINDEL>* indel_list)
{    
    FILE* fp_indel = fopen(indelfile, "r");
    if(fp_indel == NULL)                                // double check
    {
        printf("Cannot open file of indels. Exited. \n");
        exit(1);
    }
    (*indel_list).clear();
    while(!feof(fp_indel))
    {
        /* Pep2-1 Aa.chr1 5798335 5798336 2 TT core_nonrep 79 0.866667 1.0101 */
        char pronam[32];
        char chroid[32];
        char indbeg[32];
        char indend[32];
        char indlen[32];
        char indseq[1024]; // caution
        char sptype[256];  // support type: core_nonrep etc.
        char spnumb[32];   // support == coverage
        char concor[32];   
        char avghit[32];
        
        if(fscanf(fp_indel, "%s\n", pronam));
        if(fscanf(fp_indel, "%s\n", chroid));
        if(fscanf(fp_indel, "%s\n", indbeg));
        if(fscanf(fp_indel, "%s\n", indend));
        if(fscanf(fp_indel, "%s\n", indlen));// not used
        if(fscanf(fp_indel, "%s\n", indseq));
        if(fscanf(fp_indel, "%s\n", sptype));// not used
        if(fscanf(fp_indel, "%s\n", spnumb));
        if(fscanf(fp_indel, "%s\n", concor));
        if(fscanf(fp_indel, "%s\n", avghit));// not used
        
        unsigned long idb = atol(indbeg);
        unsigned long ide = atol(indend);
        /* caution: reg_chromosome, reg_begin, reg_end are global variables*/
        if((string)chroid==reg_chromosome && idb>=reg_begin && ide<=reg_end)
        {
            myINDEL tempINDEL;
            /* initialize a INDEL     */
            init_INDEL(&tempINDEL);
            /* update info of a INDEL */
            tempINDEL.ecotype    += (string)pronam;  // 0.project name
            tempINDEL.chromosome += (string)chroid;  // 1.chromosome
            tempINDEL.begin       = idb;             // 2.begin position
            tempINDEL.end         = ide;             // 3.end position
            tempINDEL.seq        += (string)indseq;  // 4.INDEL sequence
            tempINDEL.support     = atol(spnumb);    // 5.support/coverage
            tempINDEL.concordance = atof(concor);    // 6.concordance
            
            tempINDEL.peak_distance = idb - 0;       // suppose peak at position 0 - added by SHQ
            /* record the INDEL       */
            (*indel_list).insert(std::pair<unsigned long, myINDEL>(idb, tempINDEL));
        }
    }
    fclose(fp_indel);
}
/* initialization of a INDEL struct: 

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
*/
