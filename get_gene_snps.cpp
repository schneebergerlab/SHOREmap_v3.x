/* this function set up the CDS sequence for a specific gene for checking protein changes */

#include <map>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include "globals.h"

using namespace std;

/* singleGENEann   = gene_ann[gene_name]; snp_list from get_SNPlist(...) */
/* gene_coding_ann = coding_ann[gene_name] */
void get_gene_snps(std::string chromosome, 
                   std::string gene_name,
                   QUARTET singleGENEann, 
                   map<unsigned long, QUARTET> gene_coding_ann,
                   std::string isoform,
                   multimap<unsigned long, mySNP> snp_list,
                   GeneSNPlist* gene,
                   FILE* fplog)
{
    /* initialization */
    (*gene).gene_id        = gene_name;                         //1 std::string
    (*gene).isoform        = isoform;                           //  std::string
    (*gene).chromosome     = chromosome;                        //  std::string
    (*gene).start          = singleGENEann.start;               //  unsigned long
    (*gene).end            = singleGENEann.end;                 //5 unsigned long
    (*gene).cds_length     = 0;                                 //  unsigned long
    (*gene).protein_length = 0;                                 //  unsigned long
    (*gene).orientation    = singleGENEann.orien;               //  std::string
    (*gene).SNP_lists.insert(snp_list.begin(), snp_list.end()); //  map<unsigned long, mySNP>
    (*gene).coding_SNP.clear();                            //10 map<unsigned long, unsigned long>
    (*gene).NS_changes.clear();                            //   map<unsigned long, unsigned long>
    (*gene).AA_changes.clear();                            //   map<unsigned long, unsigned long>
    (*gene).CDSexon.clear();                               //   map<unsigned long, unsigned long>
    (*gene).protein.clear();                               //   map<std::string, std::string>
    (*gene).ref_gene       = singleGENEann.seq;            //15 std::string
    (*gene).ref_coding.clear();                            //   std::string
    (*gene).eco_coding.clear();                            //17 std::string
    
    /* extract and concatenate the coding region as one sequence */
    map<unsigned long, QUARTET>::iterator coding_itr;
    map<unsigned long, QUARTET>::iterator coding_itr_end;
    coding_itr     = gene_coding_ann.begin();
    coding_itr_end = gene_coding_ann.end();
    unsigned long icoding_start = 999999999;
    unsigned long icoding_end   = 0;
    // first_codon_frame
    // last_codon_frame
    while(coding_itr != coding_itr_end)
    {
        /* start and end positions for each CDS */
        std::pair <unsigned long, unsigned long> cds_pair;
        cds_pair = std::make_pair((*coding_itr).second.start, (*coding_itr).second.end);
        (*gene).CDSexon.insert(cds_pair);
        /* increase length of overall CDS */
        (*gene).cds_length += (*coding_itr).second.end - (*coding_itr).second.start + 1;
        /* concatenate all the sequences of CDS */
        (*gene).ref_coding += (*coding_itr).second.seq;
        // get the frame value of the first and last codon
        if((*coding_itr).second.start < icoding_start) 
        {
            icoding_start     = (*coding_itr).second.start;
            first_codon_frame = (*coding_itr).second.frame; 
        }
        if((*coding_itr).second.end  > icoding_end) 
        {
            icoding_end       = (*coding_itr).second.end;
            last_codon_frame  = (*coding_itr).second.frame; 
        }
        coding_itr ++;
    }
    /* reverse the above sequence if orientation is "-" */  
    if((*gene).orientation == "-")
    {
        std::string cseq = (*gene).ref_coding;
        /* reverse sequence    */
        std::reverse(cseq.begin(), cseq.end());
        /* lowercase sequence  */
        std::transform(cseq.begin(), cseq.end(), cseq.begin(), ::tolower);
        /* complement sequence */
        std::replace(cseq.begin(), cseq.end(), 'a', 'T');
        std::replace(cseq.begin(), cseq.end(), 'c', 'G'); // caution1: here c:G
        std::replace(cseq.begin(), cseq.end(), 'g', 'C'); // caution2: here g:C
        std::replace(cseq.begin(), cseq.end(), 't', 'A'); // if caution 1&2 is not consistent,
        std::replace(cseq.begin(), cseq.end(), 'u', 'A'); // they will result in c:G+G:C=c:C
        /* record the reversed sequence */
        (*gene).ref_coding.clear();
        (*gene).ref_coding = cseq;
    }
    // find the exact starting point of translation
    if((*gene).orientation == "-")
    {
        first_codon_frame  = last_codon_frame;
    }
    
    if((*gene).cds_length%3 != 0) 
    {
        /* TODO: translate with shift 1/2/3, as this does happen as seen from the following cmd */
        printf("Warning (in get_gene_snps.cpp): length of cds is not a multiple of 3: ");
        printf("%s: gene range = %ld~%ld\n", (*gene).gene_id.c_str(), (*gene).start, (*gene).end);
        printf("Phase of first codon is %d, meaning the first %d base(s) in CDS will be discarded in translation. \n", first_codon_frame, first_codon_frame);
        fprintf(fplog, "Warning (in get_gene_snps.cpp): length of cds is not a multiple of 3: ");
        fprintf(fplog, "%s: gene range = %ld~%ld\n", (*gene).gene_id.c_str(), (*gene).start, (*gene).end);
        fprintf(fplog, "Phase of first codon is %d, meaning the first %d base(s) in CDS will be discarded in translation. \n", first_codon_frame, first_codon_frame);        
    }
    (*gene).protein_length = ((*gene).cds_length-first_codon_frame)/3;
}
