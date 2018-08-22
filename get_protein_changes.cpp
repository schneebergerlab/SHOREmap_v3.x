/* this function calculates CDS&proteins changes of all SNPs in tSNPlist                           */
#include   <stdio.h>
#include       <map>
#include <algorithm>
#include    <string>

#include "globals.h"
#include "translate_as_string.h"

void get_protein_changes(GeneSNPlist* gene, multimap<unsigned long, mySNP>* tSNPlist)
{
    /* 1.translate the reference CDS sequence */
    std::string pro_seq("");
    translate_as_string((*gene).ref_coding.substr(first_codon_frame), &pro_seq);
    (*gene).protein.insert(std::pair<std::string, std::string>("ref", pro_seq));
    
    /* 2.construct the mutant CDS sequence that will be translated into protein */
    unsigned long cds_pos = 1; // ID
    map<unsigned long, unsigned long>::iterator cds_exon_itr;
    map<unsigned long, unsigned long>::iterator cds_exon_itr_end;
    cds_exon_itr     = (*gene).CDSexon.begin();
    cds_exon_itr_end = (*gene).CDSexon.end();
    while(cds_exon_itr != cds_exon_itr_end)
    {
        unsigned long icds_exon = (*cds_exon_itr).first;
        while(icds_exon <= (*cds_exon_itr).second) // continuous, e.g., first=2533 ~ second=2789
        {
            std::string ref_base = (*gene).ref_gene.substr(icds_exon - (*gene).start, 1);
            
            map<unsigned long, mySNP>::iterator posSNP_itr = (*gene).SNP_lists.find(icds_exon);
            if( posSNP_itr != (*gene).SNP_lists.end() )
            {
                (*posSNP_itr).second.gene_id = (*gene).gene_id; // (*posSNP_itr) is mySNP struct
                (*posSNP_itr).second.stype   = "CDS";
                if((*gene).orientation == "+")
                {
                    std::pair <unsigned long, unsigned long> cds_pair;
                    // first=virtual(1,2,3,...,cds_length), second=real(along genome)
                    cds_pair = std::make_pair(cds_pos, icds_exon); 
                    (*gene).coding_SNP.insert(cds_pair);
                    // virtual ID
                    (*posSNP_itr).second.cds_pos  = cds_pos;
                    // relative position - caution here: correctness not checked - 2013-04-15 23:57
                    (*posSNP_itr).second.gene_pos = icds_exon - (*gene).start + 1;
                }
                else if((*gene).orientation == "-")
                {
                    // caution here: correctness not checked - 2013-04-16 00:04
                    unsigned long real_pos = (*gene).cds_length - cds_pos + 1; 
                    std::pair <unsigned long, unsigned long> cds_pair;
                    cds_pair = std::make_pair(real_pos, icds_exon);
                    (*gene).coding_SNP.insert(cds_pair);
                    
                    (*posSNP_itr).second.cds_pos  = real_pos;
                    (*posSNP_itr).second.gene_pos = (*gene).end - icds_exon + 1;
                }
                else
                {
                    printf("Warning (get_protein_changes.cpp): orient neither '+' nor '-'. ");
                    printf("Check GFF format.\n");
                }
                (*gene).eco_coding += (*posSNP_itr).second.new_base;
            }  
            else
            {
                (*gene).eco_coding += ref_base;
            }  
            /* next position of an exon */     
            cds_pos  ++;
            icds_exon ++;
        }
        /* next exon */
        cds_exon_itr ++;
    }
    
    if((*gene).orientation == "-")
    {
        std::string cseq = (*gene).eco_coding;
        /* reverse sequence    */
        std::reverse(cseq.begin(), cseq.end());
        /* lowercase sequence  */
        std::transform(cseq.begin(), cseq.end(),cseq.begin(), ::tolower);
        /* complement sequence */
        std::replace(cseq.begin(), cseq.end(), 'a', 'T');
        std::replace(cseq.begin(), cseq.end(), 'c', 'G'); // caution1: here c:G
        std::replace(cseq.begin(), cseq.end(), 'g', 'C'); // caution2: here g:C
        std::replace(cseq.begin(), cseq.end(), 't', 'A'); // if caution 1&2 is not consistent,
        std::replace(cseq.begin(), cseq.end(), 'u', 'A'); // they will result in c:G+G:C=c:C
        
        (*gene).eco_coding.clear();
        (*gene).eco_coding = cseq;
    }
    
    /* 3.translate the mutant CDS sequence */
    std::string pro_mut_seq("");
    translate_as_string((*gene).eco_coding.substr(first_codon_frame), &pro_mut_seq);
    (*gene).protein.insert(std::pair<std::string, std::string>("alt", pro_mut_seq));
    
    /* 4.get protein changes position by position */
    for(unsigned long j_prot_pos = 1; j_prot_pos <= (*gene).protein["alt"].length(); j_prot_pos ++)
    {
        map<std::string, std::string>::iterator protein_itr;
        protein_itr = (*gene).protein.find("ref");
        std::string ref_aa = (*protein_itr).second.substr(j_prot_pos-1, 1);
        protein_itr = (*gene).protein.find("alt");
        std::string new_aa = (*protein_itr).second.substr(j_prot_pos-1, 1);
        
        map<unsigned long, unsigned long>::iterator codingSNP_itr;
        /* note: protein position <=> DNA 3 positions:  1=1/2/3, 2=4/5/6, 3=7/8/9, 4=10/11/12, ...                  */
        /* there is a possibility for these positions to be SNPs                                                    */
        /* check if a SNP position 3*j_prot_pos-2+first_codon_frame exists for protein position j_prot_pos          */
        codingSNP_itr = (*gene).coding_SNP.find(3*j_prot_pos - 2 + first_codon_frame);
        if(codingSNP_itr != (*gene).coding_SNP.end())
        {
            unsigned long genome_position = (*codingSNP_itr).second;   // real position along genome
            map<unsigned long, mySNP>::iterator gpos_itr = (*gene).SNP_lists.find(genome_position);
            if(gpos_itr == (*gene).SNP_lists.end()) 
            {
                 printf("Warning (get_protein_changes.cpp): genome position not found.\n");
            }
            else
            {
                (*gpos_itr).second.codon_pos = 1;
                (*gpos_itr).second.ref_aa    = ref_aa;
                (*gpos_itr).second.new_aa    = new_aa;
                if(ref_aa != new_aa)
                {
                    /* AA change  */
                    std::pair <unsigned long, unsigned long> chg_pair;
                    chg_pair = std::make_pair(j_prot_pos, 1);
                    (*gene).AA_changes.insert(chg_pair);
                    /* DNA change */
                    chg_pair = std::make_pair(3*j_prot_pos-2, genome_position);
                    (*gene).NS_changes.insert(chg_pair);
                    (*gpos_itr).second.ns_change = 1;
                    /* stop codon change */
                    if(new_aa == "*")
                    {
                        (*gpos_itr).second.new_stop = 1;
                    }
                    else
                    if(ref_aa == "*")
                    {
                        (*gpos_itr).second.lost_stop = 1;
                    }
                }
            } 
        } // end of SNP on 3*j_prot_pos - 2
        /* check if a SNP position 3*j_prot_pos-1 exists for protein position j_prot_pos          */
        codingSNP_itr = (*gene).coding_SNP.find(3*j_prot_pos - 1 + first_codon_frame);      //diff1
        if(codingSNP_itr != (*gene).coding_SNP.end())
        {
            unsigned long genome_position = (*codingSNP_itr).second; // real position along genome
            map<unsigned long, mySNP>::iterator gpos_itr = (*gene).SNP_lists.find(genome_position);
            if(gpos_itr == (*gene).SNP_lists.end()) 
            {
                 printf("Warning (get_protein_changes.cpp): genome position not found.\n");
            }
            else
            {
                (*gpos_itr).second.codon_pos = 2;                                           // diff2
                (*gpos_itr).second.ref_aa    = ref_aa;
                (*gpos_itr).second.new_aa    = new_aa;
                if(ref_aa != new_aa)
                {
                    std::pair <unsigned long, unsigned long> chg_pair;
                    chg_pair = std::make_pair(j_prot_pos, 1);
                    (*gene).AA_changes.insert(chg_pair);    
                    chg_pair = std::make_pair(3*j_prot_pos-1, genome_position);             // diff3
                    (*gene).NS_changes.insert(chg_pair);
                    (*gpos_itr).second.ns_change = 1;
                    if(new_aa == "*")
                    {
                        (*gpos_itr).second.new_stop = 1;
                    }
                    else
                    if(ref_aa == "*")
                    {
                        (*gpos_itr).second.lost_stop = 1;
                    }
                }
            } 
        } // end of SNP on 3*j_prot_pos - 1   
        /* check if a SNP position 3*j_prot_pos exists for protein position j_prot_pos            */
        codingSNP_itr = (*gene).coding_SNP.find(3*j_prot_pos + first_codon_frame);          //diff1
        if(codingSNP_itr != (*gene).coding_SNP.end())
        {
            unsigned long genome_position = (*codingSNP_itr).second;   // real position along genome
            map<unsigned long, mySNP>::iterator gpos_itr = (*gene).SNP_lists.find(genome_position);
            if(gpos_itr == (*gene).SNP_lists.end()) 
            {
                 printf("Warning (get_protein_changes.cpp): genome position not found.\n");
            }
            else
            {
                (*gpos_itr).second.codon_pos = 3;                                            //diff2
                (*gpos_itr).second.ref_aa    = ref_aa;
                (*gpos_itr).second.new_aa    = new_aa;
                if(ref_aa != new_aa)
                {
                    std::pair <unsigned long, unsigned long> chg_pair;
                    chg_pair = std::make_pair(j_prot_pos, 1);
                    (*gene).AA_changes.insert(chg_pair);
                    chg_pair = std::make_pair(3*j_prot_pos, genome_position);               // diff3
                    (*gene).NS_changes.insert(chg_pair);  
                    (*gpos_itr).second.ns_change = 1;
                    if(new_aa == "*")
                    {
                        (*gpos_itr).second.new_stop = 1;
                    }
                    else
                    if(ref_aa == "*")
                    {
                        (*gpos_itr).second.lost_stop = 1;
                    }
                }
            } 
        } // end of SNP on 3*j_prot_pos
    }     // end of step 4. 
    (*tSNPlist).clear(); // optimization: only need to update gene-related SNPs. - 2013-04-16 21:01
    (*tSNPlist).insert( (*gene).SNP_lists.begin(), (*gene).SNP_lists.end() );
}

/* a gene struct:
    (*gene).gene_id        = gene_name;           //1   std::string
    (*gene).isoform        = isoform;             //    std::string
    (*gene).chromosome     = chromosome;          //    std::string
    (*gene).start          = singleGENEann.start; //    unsigned long
    (*gene).end            = singleGENEann.end;   //5   unsigned long
    (*gene).cds_length     = 0;                   //    unsigned long
    (*gene).protein_length = 0;                   //    unsigned long
    (*gene).orientation    = singleGENEann.orien; //    std::string
    (*gene).SNP_lists      = snp_list;            //    map<unsigned long, mySNP>
    (*gene).coding_SNP.clear();                   //10  map<unsigned long, unsigned long>
    (*gene).NS_changes.clear();                   //    map<unsigned long, unsigned long>
    (*gene).AA_changes.clear();                   //    map<unsigned long, unsigned long>
    (*gene).CDSexon.clear();                      //    map<unsigned long, unsigned long>
    (*gene).protein.clear();                      //    map<std::string, std::string>
    (*gene).ref_gene       = singleGENEann.seq;   //15  std::string
    (*gene).ref_coding.clear();                   //    std::string
    (*gene).eco_coding.clear();                   //17  std::string
*/

/* initialization of mySNP struct* iSNP: 
    (*iSNP).ecotype        = "";
    (*iSNP).chromosome     = "";
    (*iSNP).position       = 0;
    (*iSNP).stype          = "intergenic";
    (*iSNP).gene_id        = "NA";
    (*iSNP).gene_pos       = 0;
    (*iSNP).cds_pos        = 0;
    (*iSNP).codon_pos      = 0;
    (*iSNP).ns_change      = 0;
    (*iSNP).new_stop       = 0;
    (*iSNP).lost_stop      = 0;
    (*iSNP).splicechange   = 0;
    (*iSNP).ref_base       = "";
    (*iSNP).new_base       = "";
    (*iSNP).ref_aa         = "";
    (*iSNP).new_aa         = "";
    //map<std::string, std::string>domain_change;
    (*iSNP).support        = 0;
    (*iSNP).concordance    = 0.0;
    (*iSNP).quality        = 0;
    (*iSNP).peak_distance  = 999999999;
    (*iSNP).marker_ratio   = 0.0;
    (*iSNP).source         = "IlluminaGA2";
*/
