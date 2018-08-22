/* this function reads info of SNPs within a region of a specific chromosome from a SNP file      */

#include <map>
#include <stdlib.h>
#include <stdio.h>
#include "init_SNP.h"

void get_SNPlist(char* snpfile, multimap<unsigned long, mySNP>* snp_list)
{    
    FILE* fp_snp = fopen(snpfile, "r");
    if(fp_snp == NULL)                                // double check
    {
        printf("Cannot open file of SNPs. Exited. \n");
        exit(1);
    }
    (*snp_list).clear();
    while(!feof(fp_snp))
    {
        char projname[32];
        char chrid[32];
        char position[32];                            // caution: not long type.
        char refb[32];                                // caution: char
        char mutb[32];                                // caution: char
        char qscore[32];
        char pcoverage[32];
        char pconcordance[32];
        char pavghits[32];
        
        fscanf(fp_snp, "%s\n", projname);
        fscanf(fp_snp, "%s\n", chrid);
        fscanf(fp_snp, "%s\n", position);
        fscanf(fp_snp, "%s\n", refb);
        fscanf(fp_snp, "%s\n", mutb);
        fscanf(fp_snp, "%s\n", qscore);
        fscanf(fp_snp, "%s\n", pcoverage);
        fscanf(fp_snp, "%s\n", pconcordance);
        fscanf(fp_snp, "%s\n", pavghits);
        
        unsigned long posn = atol(position);
        /* caution: reg_chromosome, reg_begin, reg_end are global variables*/
        if((string)chrid==reg_chromosome && posn>=reg_begin && posn<=reg_end)
        {
            mySNP tempSNP;
            /* initialize a SNP     */
            init_SNP(&tempSNP);
            /* update info of a SNP */
            tempSNP.ecotype    += (string)projname;   // 0.project name
            tempSNP.chromosome += (string)chrid;      // 1.chromosome
            tempSNP.position    = posn;               // 2.position
            tempSNP.ref_base    = (string)refb;       // 3.ref base - caution
            tempSNP.new_base    = (string)mutb;       // 4.new base
            tempSNP.quality     = atof(qscore);       // 5.quality score
            tempSNP.support     = atol(pcoverage);    // 6.support/coverage
            tempSNP.concordance = atof(pconcordance); // 7.concordance
            /* record the SNP       */
            (*snp_list).insert(std::pair<unsigned long, mySNP>(posn, tempSNP));
        }
    }
    fclose(fp_snp);
}
/* initialization of a SNP struct: 
    (*iSNP).ecotype       = "";
    (*iSNP).chromosome    = "";
    (*iSNP).position      = 0;
    (*iSNP).stype         = "intergenic";
    (*iSNP).gene_id       = "NA";
    (*iSNP).gene_pos      = 0;
    (*iSNP).cds_pos       = 0;
    (*iSNP).codon_pos     = 0;
    (*iSNP).ns_change     = 0;
    (*iSNP).new_stop      = 0;
    (*iSNP).lost_stop     = 0;
    (*iSNP).splicechange  = 0;
    (*iSNP).ref_base      = "";
    (*iSNP).new_base      = "";
    (*iSNP).ref_aa        = "";
    (*iSNP).new_aa        = "";
    //map<std::string, std::string>domain_change;
    (*iSNP).support       = 0;
    (*iSNP).concordance   = 0.0;
    (*iSNP).quality       = 0;
    (*iSNP).peak_distance = 999999999;
    (*iSNP).marker_ratio  = 0.0;
    (*iSNP).source        = "IlluminaGA2";
*/
