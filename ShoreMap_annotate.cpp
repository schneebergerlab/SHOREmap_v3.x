/* UPDATE: 
   Date  : 2015-OCT-09 : add function for visualizing annotations of snps along chromosome               
*/
#include    <stdlib.h>
#include     <stdio.h>
#include         <map>
#include      <math.h>
#include      <string>
#include    <string.h>
#include    <iostream>
// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>
#include "globals.h"
#include "precheck_opt.h"
#include "print_error_exit.h"
#include "read_chromosomes.h"
#include "get_SNPlist.h"
#include "get_INDELlist.h"
#include "read_GFFinfo.h"
#include "get_gene_snps.h"
#include "get_protein_changes.h"
#include "read_peaks.h"
#include "compare2strings.h"
#include "plot_chr_snp_annotations.h"

void cmd_init_annotate(int argc, char* argv[]);

void ShoreMap_annotate(int argc, char* argv[])
{
    if(argc < 10) // at least the first 5 parameters provided.
    {
        printf("\n\nUsage: SHOREmap annotate [options]\n\n");
        printf("#Mandatory:\n");
        printf("--chrsizes   STRING     File of name and size of chromosomes\n");
        printf("--snp        STRING     File of mutations to annotate\n");
        printf("--chrom      STRING     ID of Chromosome to annotate\n");
        printf("--start      INT        Start of targeted region\n");
        printf("--end        INT        End of targeted region\n");
        printf("--folder     STRING     output folder.\n\n");
        printf("#Functional annotation:\n");
        printf("--genome     STRING     FASTA-File of reference sequence (chromosomes)\n");
        printf("--gff        STRING     GFF-file of gene annotation\n\n");
        printf("#Optional:\n");
        printf("--peaks      STRING     File of peaks info (for ranking mutations) \n\n");
        printf("--ins        STRING     File of insertions\n");
        printf("--del        STRING     File of deletions\n\n");
        printf("-vis         BOOL       Trun on visualization of annotation (require --consen).\n");
        printf("--consen     STRING     File of fg-consensus; necessary if -vis is on;[NULL]\n");
        //printf("SHQ-note2: codes for indels can be reduced as functions - TODO .\n");
        exit(1);
    }
    
    /* set parameters according to cmd line */
    cmd_init_annotate(argc, argv);
    /* write annotate log */
    std::string fANNlog("");
    fANNlog = out_folder +  "Annotation.log";
    FILE* fplog = fopen(fANNlog.c_str(), "a+");
    if(fplog)
    {
        fprintf(fplog, "Inputs provided:\n\n");
        map<string, string>::iterator cmditr = CMD.begin();
        while(cmditr != CMD.end())
        {
            if((*cmditr).first.length()<8)
            {
                fprintf(fplog, "%s \t\t%s\n", (*cmditr).first.c_str(), (*cmditr).second.c_str());
            }
            else
            {
                fprintf(fplog, "%s \t%s\n", (*cmditr).first.c_str(), (*cmditr).second.c_str());
            }
            cmditr++;
        }
    }
    else
    {
        printf("Cannot open file to write log. Exited (in init_backcross).\n"); exit(1);
    }
    
    /* check coding change file - can be removed after future testing - 2013-04-16 20:30          */
    char checkfile[512];
    sprintf(checkfile, "%sref_and_eco_coding_seq_%s.txt\0", out_folder.c_str(), reg_chromosome.c_str());
    FILE* fpcheck = fopen(checkfile, "w");
    if(fpcheck == NULL)
    {
        printf("Cannot open file to write coding DNA and translated protein. Exited. \n");
        exit(1);
    }
    if(verbose) cout << "File for recording ref and eco coding sequences is ready. " << endl;
    /* end preparation of checking                                                                */
    
    /* TODO: check ref error - 2013-04-16 */
    map<std::string, map<unsigned long, int> > ref_err;
    ref_err.clear();
    /* need a function here to read reference error */
    
    /* Get interested SNPs of all genes from file such as SHOREmap_marker.bg_corrected_q25_f20_EMS*/
    if(verbose) cout << "Reading SNPs from file " << snp_file << endl;
    get_SNPlist((char*)snp_file.c_str(), &SNPlist);
    if(verbose) cout << SNPlist.size() << " SNPs within the specified chr-region have been read from file: " << snp_file <<  endl;
    multimap<unsigned long, mySNP>::iterator snpitr = SNPlist.begin();
    while(false && snpitr != SNPlist.end())
    {
        printf("%s\t",   (*snpitr).second.ecotype.c_str());
        printf("%s\t",   (*snpitr).second.chromosome.c_str());
        printf("%ld\t",  (*snpitr).second.position);
        printf("%s\t",   (*snpitr).second.ref_base.c_str());
        printf("%s\t",   (*snpitr).second.new_base.c_str());
        printf("%f\t",   (*snpitr).second.quality);
        printf("%ld\t",  (*snpitr).second.support);
        printf("%f\t\n", (*snpitr).second.concordance);
        snpitr ++;
    }
    /* Get insertion, deletion list info */
    if(ins_file.length() > 0) 
    {
        INSlist.clear();
        get_INDELlist((char*)ins_file.c_str(), &INSlist);
        if(verbose)
        {
            printf("%ld insertions have been read from file: %s\n", 
                   INSlist.size(), ins_file.c_str());
        }
    }
    if(del_file.length() > 0) 
    {
        DELlist.clear();
        get_INDELlist((char*)del_file.c_str(), &DELlist);
        if(verbose)
        {
            printf("%ld deletions have been read from file: %s\n",
                   DELlist.size(), del_file.c_str());
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /* prioritize snps, indels with the distance from the peak position - TODO 2013-04-25 14:41   */
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /* Read gene annotation from file such as A_alpina_V3.annotation.gff3                         */
    std::string isoform = "";
    read_GFFinfo( (char*)gff_file.c_str(), &isoform);
    if(verbose) printf("Annotation info on %ld genes is ready.\n", gene_ann.size());
    if(gene_ann.size() == 0)
    {
        printf("\n\nERROR: for specified chromosome %s: \n", (char*)reg_chromosome.c_str());
        printf("   no gene annotations are found in %s; \n", (char*)gff_file.c_str());
        printf("   hint1: check content/format of gff file provided; \n");
        printf("   hint2: check consistency of chromosome-IDs in all files and the one to annotate.\n");
        exit(1);
    }
    /* check if there are seq_type ann exceeding given chrsizes (and display seq_type on screen)  */
    multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
    ann_itr = seq_type.begin();
    while(ann_itr != seq_type.end())
    {
        //printf("%s:\n", (*ann_itr).first.c_str());
        multimap<std::string, RANGE> typeMap = (*ann_itr).second;
        multimap<std::string, RANGE>::iterator type_itr;
        type_itr = typeMap.begin();
        while(type_itr != typeMap.end())
        {
            /* check if the length seq type is within chromosome region */
            if((*type_itr).second.start <=0 || (*type_itr).second.end > CHR2SIZE[reg_chromosome])
            {
                fprintf(fplog, "Warning: region of %s exceeds chromosome size given.\n", 
                       (char*)(*type_itr).first.c_str());
                printf("start=%d; end=%d, while chr size is %d", (*type_itr).second.start, (*type_itr).second.end, CHR2SIZE[reg_chromosome]);
                printf("Warning: %s exceeds chromosome size given. ", 
                       (char*)(*type_itr).first.c_str());
                printf("Please check the gff file and chromosome sizes file. Exited.\n");
                exit(1);
            }
            /*
            printf("\t%s\t%ld\t%ld\t\n",
                 (*type_itr).first.c_str(),
                 (*type_itr).second.start,
                 (*type_itr).second.end);
            */     
            type_itr ++;
        }
        ann_itr ++;
    }
    /*   check if there are gene_ann exceeding given chrsizes (and display gene_ann on screen)    */
    map<std::string, QUARTET>::iterator gitr;
    gitr = gene_ann.begin();
    while(gitr != gene_ann.end())
    {
        /* check if the length seq type is within chromosome region */
        if((*gitr).second.start <=0 || (*gitr).second.end > CHR2SIZE[reg_chromosome])
        {
             fprintf(fplog, "Warning: region of %s exceeds chromosome size given.\n", 
                   (char*)(*gitr).first.c_str());
             printf("Warning: %s exceeds chromosome size given. ", 
                   (*gitr).first.c_str());
             printf("Please check the gff file and chromosome sizes file. Exited.\n");
             exit(1);
        }
        /*
        printf("%s\t%ld\t%ld\t%ld\t%s\n", 
            (*gitr).first.c_str(), 
            (*gitr).second.start, 
            (*gitr).second.end, 
            (*gitr).second.seq.length(), 
            (*gitr).second.orien.c_str());
        */
        gitr ++;
    }
    /*  check if there are coding_ann exceeding given chrsizes (and display coding_ann on screen) */
    map<std::string, map<unsigned long, QUARTET> >::iterator cditr;
    cditr = coding_ann.begin();
    while(cditr != coding_ann.end())
    {
        //printf("%s:\n", (*cditr).first.c_str()); 
        map<unsigned long, QUARTET> cdMap = (*cditr).second;
        map<unsigned long, QUARTET>::iterator coitr;
        coitr = cdMap.begin();
        while(coitr != cdMap.end())
        {
            /* check if the length seq type is within chromosome region */
            if((*coitr).second.start <=0 || (*coitr).second.end > CHR2SIZE[reg_chromosome])
            {
                fprintf(fplog, "Warning: region of coding exceeds chromosome size given.\n");
                printf("Warning: coding region exceeds chromosome size given. ");
                printf("Please check the gff file and chromosome sizes file. Exited.\n");
                exit(1);
            }
            /*
            printf("\t%ld\t%ld\t%ld\t%s\n",
               (*coitr).second.start, 
               (*coitr).second.end, 
               (*coitr).second.seq.length(),
               (*coitr).second.orien.c_str());
            */
            coitr ++;
        }
        cditr ++;
    }
    
    /* Add gene SNP annotation. STRUCTURES:  */
    /* 
        map<std::string, QUARTET> gene_ann;
        map<std::string, map<unsigned long, QUARTET> > coding_ann;
        multimap<std::string, multimap<std::string, RANGE> > seq_type;
        multimap<unsigned long, mySNP> SNPlist;
    */
    
    /* SNPs should be given isoform info: caution here */
    char gene_snp_file[1024];
    sprintf(gene_snp_file, "%sgene_snp_check.txt\0", out_folder.c_str());
    FILE* fpgenesnp; 
    
    multimap<unsigned long, mySNP> SNPlist_tmp; // one position -> many SNPs due to isoforms
    multimap<unsigned long, myINDEL>    INSlist_tmp;
    multimap<unsigned long, myINDEL>    DELlist_tmp;
    SNPlist_tmp.clear();
    INSlist_tmp.clear();
    DELlist_tmp.clear();
    map<std::string, QUARTET>::iterator gene_itr;
    map<std::string, QUARTET>::iterator gene_itr_end;
    gene_itr     = gene_ann.begin();
    gene_itr_end = gene_ann.end();
    while(gene_itr != gene_itr_end) /* (genes with many isoforms are deemed as different ones)    */
    {
        /* TODO: reduce code for getting SNPs&INDELS of a gene -> function                        */
        /* SNPs of current gene - new */ 
        multimap<unsigned long, mySNP> SNPlist_curGene;
        multimap<unsigned long, mySNP>::iterator snp_itr_tot;     
        multimap<unsigned long, mySNP>::iterator snp_itr_tot_end;
        snp_itr_tot     = SNPlist.begin();// iterate SNPlist to find snps within current gene region
        snp_itr_tot_end = SNPlist.end();
        while(snp_itr_tot != snp_itr_tot_end)
        {
            /* record a SNP if it is located in this gene   */
            if((*snp_itr_tot).first>=(*gene_itr).second.start)
            if((*snp_itr_tot).first<=(*gene_itr).second.end)
            {
                std::pair<unsigned long, mySNP> snp_pair;
                snp_pair = std::make_pair((*snp_itr_tot).first, (*snp_itr_tot).second);
                /* pre-set basic gene info to a snp      */
                snp_pair.second.gene_id = (*gene_itr).first;
                snp_pair.second.stype   = "intronic/noncoding";     // can be changed to "CDS", etc.
                SNPlist_curGene.insert(snp_pair);
            }
            snp_itr_tot ++;
        }
        /* INSERTION of current gene */
        multimap<unsigned long, myINDEL> INSlist_curGene;
        multimap<unsigned long, myINDEL>::iterator ins_itr;
        multimap<unsigned long, myINDEL>::iterator ins_itr_end;
        ins_itr     = INSlist.begin();
        ins_itr_end = INSlist.end();
        while(ins_itr != ins_itr_end)
        {
            if((*ins_itr).second.begin>=(*gene_itr).second.start)
            if((*ins_itr).second.end  <=(*gene_itr).second.end)
            {
                std::pair<unsigned long, myINDEL> ins_pair;
                ins_pair = std::make_pair((*ins_itr).first, (*ins_itr).second);
                ins_pair.second.gene_id = (*gene_itr).first;
                ins_pair.second.stype   = "intronic/noncoding";
                INSlist_curGene.insert(ins_pair);
            }
            ins_itr ++;
        }
        /* DELETION of current gene */
        multimap<unsigned long, myINDEL> DELlist_curGene;
        multimap<unsigned long, myINDEL>::iterator del_itr;
        multimap<unsigned long, myINDEL>::iterator del_itr_end;
        del_itr     = DELlist.begin();
        del_itr_end = DELlist.end();
        while(del_itr != del_itr_end)
        {
            if((*del_itr).second.begin>=(*gene_itr).second.start)
            if((*del_itr).second.end  <=(*gene_itr).second.end)
            {
                std::pair<unsigned long, myINDEL> del_pair;
                del_pair = std::make_pair((*del_itr).first, (*del_itr).second);
                del_pair.second.gene_id = (*gene_itr).first;
                del_pair.second.stype   = "intronic/noncoding";
                DELlist_curGene.insert(del_pair);
            }
            del_itr ++;
        }
        
        /* Why there is no check of indel within coding region? TODO 2013-05-08 20:29?            */
        /* coding region: get the corresponding coding annotation for this gene                   */
        map<std::string, map<unsigned long, QUARTET> >::iterator coding_itr;
        map<std::string, map<unsigned long, QUARTET> >::iterator coding_itr_end;
        coding_itr     = coding_ann.find((*gene_itr).first);                        // isoform added
        coding_itr_end = coding_ann.end();
        if(coding_itr != coding_itr_end)
        {
            /* s1.get interested SNPs in (global) [region_begin, reg_end] of reg_chromosome       */
            unsigned long isofpos = (*gene_itr).first.find(".")+1;
            unsigned long isolen  = (*gene_itr).first.length() - 1 - isofpos + 1;
            std::string isoform   = (*gene_itr).first.substr(isofpos, isolen);  // caution: standard
            struct GeneSNPlist gene;
            get_gene_snps(reg_chromosome, (*gene_itr).first, (*gene_itr).second, 
                   (*coding_itr).second, isoform, SNPlist_curGene, &gene, fplog);            //  new                  
            /* s2.get protein changes */
            get_protein_changes(&gene, &SNPlist_curGene); // SNPs of this gene               //  new       

            /* check */
            fprintf(fpcheck, ">%s_orien=\'%s\'\n", gene.gene_id.c_str(), gene.orientation.c_str());
      
            string diff_str;
            compare2strings(gene.ref_coding.substr(first_codon_frame), gene.eco_coding.substr(first_codon_frame), &diff_str);
            fprintf(fpcheck, "%s\n%s\n", gene.ref_coding.substr(first_codon_frame).c_str(), diff_str.c_str());

            compare2strings(gene.protein["ref"], gene.protein["alt"].c_str(), &diff_str);
            fprintf(fpcheck, "%s\n%s\n", gene.protein["ref"].c_str(), diff_str.c_str());
        }
        
        /* non-coding region */
        multimap<std::string, multimap<std::string, RANGE> >::iterator stype_itr;
        multimap<std::string, multimap<std::string, RANGE> >::iterator stype_itr_end;
        stype_itr        = seq_type.find((*gene_itr).first);
        stype_itr_end    = seq_type.end();
        if(stype_itr != stype_itr_end)           // seq_type = <gene_name, <"CDS/gene/...", RANGE> >
        {
            /* (*stype_itr).second == multimap<"CDS/gene/...", RANGE>        <= range of seq_type */
            multimap<std::string, RANGE>::iterator type2_itr;
            multimap<std::string, RANGE>::iterator type2_itr_end;
            type2_itr_end = (*stype_itr).second.end();
            for(type2_itr = (*stype_itr).second.begin(); type2_itr != type2_itr_end; type2_itr ++)
            {
                /* (*type2_itr).first:  sequence type = "CDS/gene/exon/mRNA/..."                  */
                /* (*type2_itr).second: range of sequence type                                    */
                /* range of CDS/gene/exon/mRNA/...                                                */
                unsigned long type2_reg_start = (*type2_itr).second.start;
                unsigned long type2_reg_end   = (*type2_itr).second.end;
                
                /* SNP */
                multimap<unsigned long, mySNP>::iterator snp_itr;
                multimap<unsigned long, mySNP>::iterator snp_itr_end;
                snp_itr_end = SNPlist_curGene.end();
                for(snp_itr = SNPlist_curGene.begin(); snp_itr != snp_itr_end; snp_itr ++)    // new           
                {
                    unsigned long mysnppos = (*snp_itr).first;
                    multimap<unsigned long, mySNP>::iterator titr_snp;
                    titr_snp = SNPlist_curGene.find(mysnppos);
                    if(mysnppos>=type2_reg_start)                   // caution: overlapped gene anns
                    if(mysnppos<=type2_reg_end)
                    if((*type2_itr).first == "CDS")
                    {
                        (*titr_snp).second.stype = "CDS";
                    }
                    else if((*type2_itr).first.substr(0, 18) == "splice_site_change")     // caution
                    {
                        if((*titr_snp).second.stype != "CDS")
                        {
                            (*titr_snp).second.stype   = "splice_site_change";
                            (*titr_snp).second.gene_id = (*gene_itr).first;             // gene name
                        }
                    }
                    else if((*type2_itr).first=="five_prime_UTR"||(*type2_itr).first=="three_prime_UTR")
                    {
                        if((*titr_snp).second.stype != "CDS")
                        if((*titr_snp).second.stype != "splice_site_change")
                        {
                            (*titr_snp).second.stype   = (*type2_itr).first;
                            (*titr_snp).second.gene_id = (*gene_itr).first;             // gene name
                        }
                    }
                }
                /* INSERTION TODO require modification on isoform info 2013-05-08 18:39           */
                multimap<unsigned long, myINDEL>::iterator ins_itr;
                if(ins_file.length()>0 && INSlist_curGene.size()>0)
                for(ins_itr = INSlist_curGene.begin(); ins_itr != INSlist_curGene.end(); ins_itr ++)
                {
                   unsigned long myinsbeg = (*ins_itr).second.begin;
                   unsigned long myinsend = (*ins_itr).second.end;
                   
                   if( (myinsbeg>=type2_reg_start && myinsbeg<=type2_reg_end) || 
                       (myinsend>=type2_reg_start && myinsend<=type2_reg_end) )
                   if((*type2_itr).first == "CDS")
                    {
                        (*ins_itr).second.stype = "CDS";
                    }
                    else if((*type2_itr).first.substr(0, 18) == "splice_site_change")     // caution
                    {
                        if((*ins_itr).second.stype != "CDS")
                        {
                            (*ins_itr).second.stype   = "splice_site_change";
                            (*ins_itr).second.gene_id = (*gene_itr).first;              // gene name
                        }
                    }
                    else if((*type2_itr).first=="five_prime_UTR"||(*type2_itr).first=="three_prime_UTR")
                    {
                        if((*ins_itr).second.stype != "CDS")
                        if((*ins_itr).second.stype != "splice_site_change")// caution: diff from origin
                        {
                            (*ins_itr).second.stype   = (*type2_itr).first;
                            (*ins_itr).second.gene_id = (*gene_itr).first;              // gene name
                        }
                    }
                }
                /* DELETION TODO require modification on isoform info 2013-05-08 18:39            */
                multimap<unsigned long, myINDEL>::iterator del_itr;
                if(del_file.length()>0 && DELlist_curGene.size()>0)
                for(del_itr = DELlist_curGene.begin(); del_itr != DELlist_curGene.end(); del_itr ++)
                {
                   unsigned long mydelbeg = (*del_itr).second.begin;
                   unsigned long mydelend = (*del_itr).second.end;
                   
                   if( (mydelbeg>=type2_reg_start && mydelbeg<=type2_reg_end) ||
                       (mydelend>=type2_reg_start && mydelend<=type2_reg_end) ||
                       (mydelbeg<type2_reg_start && type2_reg_end<mydelend) )
                   if((*type2_itr).first == "CDS")
                    {
                        (*del_itr).second.stype = "CDS";
                    }
                    else if((*type2_itr).first.substr(0, 18) == "splice_site_change")     // caution
                    {
                        if((*del_itr).second.stype != "CDS")
                        {
                            (*del_itr).second.stype   = "splice_site_change";
                            (*del_itr).second.gene_id = (*gene_itr).first;              // gene name
                        }
                    }
                    else if((*type2_itr).first=="five_prime_UTR"||(*type2_itr).first=="three_prime_UTR")
                    {
                        if((*del_itr).second.stype != "CDS")
                        if((*del_itr).second.stype != "splice_site_change")// caution: diff from origin
                        {
                            (*del_itr).second.stype   = (*type2_itr).first;
                            (*del_itr).second.gene_id = (*gene_itr).first;              // gene name
                        }
                    }
                }
            }
        }  
        // new
        SNPlist_tmp.insert(SNPlist_curGene.begin(), SNPlist_curGene.end());
        INSlist_tmp.insert(INSlist_curGene.begin(), INSlist_curGene.end());
        DELlist_tmp.insert(DELlist_curGene.begin(), DELlist_curGene.end());        
        /* next gene (annotation) */
        gene_itr ++;
    }
    /* check */
    fclose(fpcheck);
    
    // new 2013-05-07 14:32
    multimap<unsigned long, mySNP>::iterator snp_itr;
    snp_itr = SNPlist.begin();
    while(snp_itr != SNPlist.end())                                // remove snp info in gene region
    {
        if(SNPlist_tmp.find((*snp_itr).first) == SNPlist_tmp.end())
        {
            snp_itr ++;
        }
        else
        {
            multimap<unsigned long, mySNP>::iterator snp_itr_tmp = snp_itr;
            snp_itr ++;
            SNPlist.erase(snp_itr_tmp);
        }
    }
    SNPlist.insert(SNPlist_tmp.begin(), SNPlist_tmp.end());   // add updated snp info in gene region
    
    /* read in peaks if provided */
    multimap<std::string, unsigned long> mpeaks;
    if(fpeaks.size() > 0)
    {
        read_peaks((char*)fpeaks.c_str(), &mpeaks);
    }
    else
    {
        mpeaks.insert(std::pair<std::string, unsigned long>(reg_chromosome, 1));
    }
    multimap<std::string, unsigned long>::iterator peak_itr;
    multimap<std::string, unsigned long>::iterator peak_itr_end;
    peak_itr     = mpeaks.begin();
    peak_itr_end = mpeaks.end();
    while(peak_itr != peak_itr_end)
    {
        /* skip peaks not located in the provided chromosome                                      */
        if((*peak_itr).first != reg_chromosome) 
        {
            peak_itr ++;
            continue; 
        }
        
        unsigned long tmp_peak = (*peak_itr).second;
        
        SNPlist_tmp.clear();
        multimap<unsigned long, mySNP>::iterator tmp_snp_itr;
        multimap<unsigned long, mySNP>::iterator tmp_snp_itr_end;
        tmp_snp_itr     = SNPlist.begin();
        tmp_snp_itr_end = SNPlist.end();
        while(tmp_snp_itr != tmp_snp_itr_end)
        {
            if((*tmp_snp_itr).first >= tmp_peak) 
            {
                unsigned long rank_pos = (*tmp_snp_itr).first - tmp_peak;
                SNPlist_tmp.insert(std::pair<unsigned long, mySNP>(rank_pos, (*tmp_snp_itr).second));
            }
            else
            {
                unsigned long rank_pos = tmp_peak - (*tmp_snp_itr).first;
                SNPlist_tmp.insert(std::pair<unsigned long, mySNP>(rank_pos, (*tmp_snp_itr).second));
            }
            tmp_snp_itr ++;
        }
        
        /* output annotations on SNPs to file */
        char out_snp_ann[1024];
        sprintf(out_snp_ann, "%sprioritized_snp_%s_%ld_%ld_peak%ld.txt\0", 
            out_folder.c_str(), reg_chromosome.c_str(), reg_begin, reg_end, tmp_peak);
        FILE* fpSNP = fopen(out_snp_ann, "w");
        if(fpSNP == NULL) 
        {
            printf("Cannot open file to write SNP annotations. Exited.\n");
            exit(1);
        }
        // new end
        // for(snp_itr = SNPlist.begin(); snp_itr != SNPlist.end(); snp_itr ++)
        for(snp_itr = SNPlist_tmp.begin(); snp_itr != SNPlist_tmp.end(); snp_itr ++)
        {
            fprintf(fpSNP, "%s\t",   reg_chromosome.c_str());                         // 1.chr
            fprintf(fpSNP, "%ld\t",  (*snp_itr).second.position);                     // 2.position
            fprintf(fpSNP, "%s\t",   (*snp_itr).second.ref_base.c_str());             // 3.ref base
            fprintf(fpSNP, "%s\t",   (*snp_itr).second.new_base.c_str());             // 4.new base
            fprintf(fpSNP, "%ld\t",  (*snp_itr).second.support);                      // 5.support
            fprintf(fpSNP, "%.2f\t", (*snp_itr).second.concordance);                  // 6.concordance
            fprintf(fpSNP, "%.0f",   (*snp_itr).second.quality);                      // 7.quality
        
            /* caution: reading of ref error has not been finished. - 2013-04-16 */
            map<std::string, map<unsigned long, int> >::iterator ireg_chr;
            ireg_chr = ref_err.find(reg_chromosome);
            if(ireg_chr!=ref_err.end()
               && (*ireg_chr).second.find( (*snp_itr).first ) != (*ireg_chr).second.end())
            {
                fprintf(fpSNP, "\tREFERR");                                           // 8.
            }
            else
            {
                fprintf(fpSNP, "\tNEWSNP");                                           // 8.
            }
        
            if(gff_file!="" && genome_file!="")
            {
                /* : "CDS", "intronic/noncoding", "intergenic"                                    */
                fprintf(fpSNP, "\t%s", (*snp_itr).second.stype.c_str());  // 9.region of SNP
                if((*snp_itr).second.gene_id != "NA")                     // 10.name of gene.isoform
                {
                    fprintf(fpSNP, "\t%s\t", (*snp_itr).second.gene_id.c_str());     
                    fprintf(fpSNP, "%ld",  (*snp_itr).first);                       // 11. dis2peak
                }
                if((*snp_itr).second.cds_pos != 0)
                {
                    std::string mySyn = "Syn";
                    if((*snp_itr).second.ns_change == 1) mySyn = "Nonsyn";
                    fprintf(fpSNP, "\t%ld", (*snp_itr).second.cds_pos);               // 12.
                    fprintf(fpSNP, "\t%ld", (*snp_itr).second.codon_pos);             // 13.
                    fprintf(fpSNP, "\t%s",  mySyn.c_str());                           // 14.
                    fprintf(fpSNP, "\t%s",  (*snp_itr).second.ref_aa.c_str());        // 15.
                    fprintf(fpSNP, "\t%s",  (*snp_itr).second.new_aa.c_str());        // 16.
                }
            }
            fprintf(fpSNP, "\n");
        }
        fclose(fpSNP);
        fprintf(fplog, "\nAnnotation of SNPs done: %s.\n", out_snp_ann);
        peak_itr ++;
    }
    
    /* output insertion annotation */
    if(ins_file.length() > 0)
    {
        // new 2013-05-08 20:34
        multimap<unsigned long, myINDEL>::iterator ins_itr;
        ins_itr = INSlist.begin();
        while(ins_itr != INSlist.end())                            // remove ins info in gene region
        {
            if(INSlist_tmp.find((*ins_itr).first) == INSlist_tmp.end())
            {
                ins_itr ++;
            }
            else
            {
                multimap<unsigned long, myINDEL>::iterator ins_itr_tmp = ins_itr;
                ins_itr ++;
                INSlist.erase(ins_itr_tmp);
            }
        }
        INSlist.insert(INSlist_tmp.begin(), INSlist_tmp.end());// add updated ins info in gene region
        // new end
        char out_ins_ann[512];
        sprintf(out_ins_ann, "%sinsertion_%s_%ld_%ld.txt\0", 
                out_folder.c_str(), reg_chromosome.c_str(), reg_begin, reg_end);
        FILE* fpINS = fopen(out_ins_ann, "w");
        if(fpINS == NULL) 
        {
            printf("Cannot open file to write annotations on insertion. Exited.\n");
            exit(1);
        }
        for(ins_itr = INSlist.begin(); ins_itr != INSlist.end(); ins_itr ++)
        {
            fprintf(fpINS, "%s\t",   reg_chromosome.c_str());              // 1.chr
            fprintf(fpINS, "%ld\t",  (*ins_itr).second.begin);             // 2.begin
            fprintf(fpINS, "%ld\t",  (*ins_itr).second.end);               // 3.end
            fprintf(fpINS, "%s\t",   (*ins_itr).second.seq.c_str());       // 4.seq
            fprintf(fpINS, "%ld\t",  (*ins_itr).second.peak_distance);     // 5.distance from peak
            fprintf(fpINS, "%ld\t",  (*ins_itr).second.support);           // 6.support
            fprintf(fpINS, "%.2f\t", (*ins_itr).second.concordance);       // 7.concordance
            fprintf(fpINS, "%.0f",   (*ins_itr).second.quality);           // 8.quality 
            if(gff_file.length()>0 && genome_file.length()>0)
            {
                fprintf(fpINS, "\t%s", (*ins_itr).second.stype.c_str());
                if((*ins_itr).second.gene_id != "NA")
                {
                    fprintf(fpINS, "\t%s\t1", (*ins_itr).second.gene_id.c_str());
                }
            }
            fprintf(fpINS, "\n");                       
        }
        fclose(fpINS);
        fprintf(fplog, "Annotation of insertions done: %s.\n", out_ins_ann);
    }   
    /* output deletion annotation */
    if(del_file.length() > 0)
    {
        // new 2013-05-08 21:03
        multimap<unsigned long, myINDEL>::iterator del_itr;
        del_itr = DELlist.begin();
        while(del_itr != DELlist.end())                            // remove del info in gene region
        {
            if(DELlist_tmp.find((*del_itr).first) == DELlist_tmp.end())
            {
                del_itr ++;
            }
            else
            {
                multimap<unsigned long, myINDEL>::iterator del_itr_tmp = del_itr;
                del_itr ++;
                DELlist.erase(del_itr_tmp);
            }
        }
        DELlist.insert(DELlist_tmp.begin(), DELlist_tmp.end());// add updated ins info in gene region
        // new end
        char out_del_ann[512];
        sprintf(out_del_ann, "%sdeletion_%s_%ld_%ld.txt\0", 
                out_folder.c_str(), reg_chromosome.c_str(), reg_begin, reg_end);
        FILE* fpDEL = fopen(out_del_ann, "w");
        if(fpDEL == NULL) 
        {
            printf("Cannot open file to write annotations on deletion. Exited.\n");
            exit(1);
        }
        for(del_itr = DELlist.begin(); del_itr != DELlist.end(); del_itr ++)
        {
            fprintf(fpDEL,   "%s\t", reg_chromosome.c_str());              // 1.chr
            fprintf(fpDEL,  "%ld\t", (*del_itr).second.begin);             // 2.begin
            fprintf(fpDEL,  "%ld\t", (*del_itr).second.end);               // 3.end
            fprintf(fpDEL,   "%s\t", (*del_itr).second.seq.c_str());       // 4.seq
            fprintf(fpDEL,  "%ld\t", (*del_itr).second.peak_distance);     // 5.distance from peak
            fprintf(fpDEL,  "%ld\t", (*del_itr).second.support);           // 6.support
            fprintf(fpDEL, "%.2f\t", (*del_itr).second.concordance);       // 7.concordance
            fprintf(fpDEL, "%.0f",   (*del_itr).second.quality);           // 8.quality 
            if(gff_file.length()>0 && genome_file.length()>0)
            {
                fprintf(fpDEL, "\t%s", (*del_itr).second.stype.c_str());
                if((*del_itr).second.gene_id != "NA")
                {
                    fprintf(fpDEL, "\t%s\t1", (*del_itr).second.gene_id.c_str());
                }
            }
            fprintf(fpDEL, "\n");                       
        }
        fclose(fpDEL);
        fprintf(fplog, "Annotation of deletions done: %s.\n", out_del_ann);
    }
    fprintf(fplog, "\n\nOutputs achieved.\n");
    
    // visualization of annotaiton of snps
    if(vis_annotation==1)
    {
        if(verbose) printf("Visualizing annotation of snps...\n");
        if(plot_chr_snp_annotations())
        {
            fprintf(fplog, "\nAnnotations visualized in .pdf.\n");
            if(verbose) 
            printf("Annotations visualized in .pdf.\n");
        }
        else
        {
            fprintf(fplog, "\nAnnotations asked to be visualized, however failed. Please check input files and redo annotate.\n");
            if(verbose)
            printf("Annotations asked to be visualized, however failed. Please check input files and redo annotate.\n");
        }
    }
    
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    fprintf(fplog, "\nAnnotate function successfully finished on %s\n", asctime(tinfo));
    fclose(fplog);
}

/* initialize cmdline parameters */
void cmd_init_annotate(int argc, char* argv[])
{
   /* 
      Read practical values from cmd options; option-values are checked and updated
      in map<string, string> CMD with format: CMD["option"] = "arg" (if exist).
    */
    int ic = 2; // check verbose first.
    while (ic < argc)
    {
        if(!strcmp(argv[ic],"-verbose"))             // option 0 to variable: verbose
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            verbose         = 1;
            printf("Be talkative during process. \n");
        }
        ic ++;
    } 
    ic = 2;      // option ic=0: SHOREmap; ic=1: outcross/backcross/extract/annotate...
    while (ic < argc) 
    {
        if(ic == 2 && strcmp(argv[ic],"--chrsizes"))
        {
            printf("--chrsizes must be input before other options. Exited.\n");
            exit(1);
        }
        if(!strcmp(argv[ic],"--chrsizes"))           // option 1 to variable: fchrsizes
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fchrsizes        += argv[ic];
            FILE* fp_chrsizes = fopen((char*)fchrsizes.c_str(), "r");
            if(!fp_chrsizes)  { print_error_exit(argv[ic-1], false); }
            fclose(fp_chrsizes);
            if (verbose) printf("Chrs sizes read from file:\t\t%s\n", fchrsizes.c_str());
            if(!read_chromosomes((char*)fchrsizes.c_str()))
            {
                printf("ERROR: invalid content in file %s.\n", fchrsizes.c_str());
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--snp"))           // option 2 to variable: snp_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            snp_file    += argv[ic];
            FILE* fp_snp = fopen((char*)snp_file.c_str(), "r"); 
            if(fp_snp  == NULL)
            {
                printf("SNP file \'%s\' does NOT exist. Exited.\n", snp_file.c_str());
                exit(1);
            }
            if (verbose) printf("File of SNPs provided:\t\t%s\n", snp_file.c_str());
            fclose(fp_snp);
        } // --ins
        else if(!strcmp(argv[ic],"--ins"))           // option 2.1 to variable: ins_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            ins_file    += argv[ic];
            FILE* fp_ins = fopen((char*)ins_file.c_str(), "r"); 
            if(fp_ins  == NULL)
            {
                printf("INSERTION file \'%s\' does NOT exist. Exited.\n", ins_file.c_str());
                exit(1);
            }
            if (verbose) printf("File of INSERTIONs provided:\t\t%s\n", ins_file.c_str());
            fclose(fp_ins);
        }
        else if(!strcmp(argv[ic],"--del"))           // option 2.2 to variable: del_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            del_file    += argv[ic];
            FILE* fp_del = fopen((char*)del_file.c_str(), "r"); 
            if(fp_del  == NULL)
            {
                printf("DELETION file \'%s\' does NOT exist. Exited.\n", del_file.c_str());
                exit(1);
            }
            if (verbose) printf("File of DELETIONS provided:\t\t%s\n", del_file.c_str());
            fclose(fp_del);
        }
        else if(!strcmp(argv[ic], "--chrom"))       // option 3 to variable: reg_chromosome
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            reg_chromosome     += argv[ic];
            if (verbose) printf("reg_chromosome id set as:\t\t%s.\n", reg_chromosome.c_str());
            if(CHR2SIZE.find(reg_chromosome)==CHR2SIZE.end()) 
            {
                printf("ERROR: chromosome %s to zoom is not found in the list. Exited.\n", 
                    (char*)reg_chromosome.c_str()); 
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--start"))        // option 4 to variable: reg_begin
        {
            if(reg_chromosome == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_begin      = atol(argv[ic]);
            string schr(reg_chromosome);
            if(reg_begin<=0 || reg_begin>CHR2SIZE[schr])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                     argv[ic-1]); exit(1);}
            if (verbose) printf("reg_begin set as:\t\t\t%ld\n", reg_begin);
        }
        else if(!strcmp(argv[ic],"--end"))          // option 5 to variable: reg_end
        {
            if(reg_chromosome == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_end      = atol(argv[ic]);
            if(reg_end<=0 || reg_end>CHR2SIZE[reg_chromosome])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("reg_end set as:\t\t\t\t%ld\n", reg_end);
        }
        else if(!strcmp(argv[ic],"--gff"))         // option 6 to variable: gff_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            gff_file    += argv[ic];
            FILE* fp_gff = fopen((char*)gff_file.c_str(), "r"); 
            if(fp_gff   == NULL)
            {
                printf("gff file \'%s\' does NOT exist. Exited.\n", gff_file.c_str());
                exit(1);
            }
            fclose(fp_gff);
        }
        else if(!strcmp(argv[ic],"--genome"))     // option 7 to variable: genome_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            genome_file    += argv[ic];
            FILE* fp_genome = fopen((char*)genome_file.c_str(), "r"); 
            if(fp_genome   == NULL)
            {
                printf("genome file \'%s\' does NOT exist. Exited.\n", genome_file.c_str());
                exit(1);
            }
            fclose(fp_genome);
        }
        else if(!strcmp(argv[ic],"--peaks")) // option 7.5 to variable: fpeaks
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fpeaks    += argv[ic];
            FILE* fp_peaks = fopen((char*)fpeaks.c_str(), "r"); 
            if(fp_peaks   == NULL)
            {
                printf("Peaks file \'%s\' does NOT exist. Exited.\n", fpeaks.c_str());
                exit(1);
            }
            fclose(fp_peaks);
        }
        else if(!strcmp(argv[ic],"--referr"))    // option 8 to variable: freferror
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            freferror         += argv[ic];
            FILE* fp_referror  = fopen((char*)freferror.c_str(), "r"); 
            if(fp_referror == NULL)
            { 
                printf("Ref-error file \'%s\' does NOT exist. Exited.\n", freferror.c_str()); 
                exit(1);
            }
            fclose(fp_referror);
        }
        else if(!strcmp(argv[ic],"--folder"))   // option 2 to variable: out_folder
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            out_folder         += argv[ic];
            if(out_folder[out_folder.length()-1] != '/') out_folder += "/";
            DIR* dir_out_folder = opendir((char*)out_folder.c_str());
            if(dir_out_folder  == NULL)
            {
                if(!mkdir((char*)out_folder.c_str(), S_IRWXU|S_IRWXG|S_IRWXO))
                {
                 /* if !mkdir() is TRUE: a new directory has to be created.
                    if out_folder is "outfolder/" instead of "/your/output/path/outfolder", 
                    an "outfolder/" will be created under the current working directory */
                    if (verbose) printf("Folder\t\t\t\t\t%s created.\n", out_folder.c_str());
                }
                else 
                {
                    printf("ERROR: cannot create output path. Exited.\n");
                    exit(1);
                }
            }
            closedir(dir_out_folder);
        }
        else if(!strcmp(argv[ic],"-vis"))                     // option to variable: vis_annotation
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            vis_annotation  = 1;
            printf("Visualization of annotation is on. \n");
        }
        else if(!strcmp(argv[ic],"--consen"))                 // option to variable: fconsensus
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fconsensus         = argv[ic];
            FILE* fp_consensus = fopen((char*)fconsensus.c_str(), "r"); 
            if(fp_consensus   == NULL)
            {
                printf("consensus file \'%s\' does NOT exist. Exited.\n", fconsensus.c_str());
                exit(1);
            }
            fclose(fp_consensus);
        }
        else if(!strcmp(argv[ic],"-verbose")) ;
        else 
        {
            printf(" cmd %s is not recognized in ShoreMap_annotate(...).\n", argv[ic]);
        }
        ic ++;
    }
    
    /* check if are necessary inputs are provided. */
    map<std::string, std::string>::iterator cmd_itr_end = CMD.end(); 
    bool cmd_missing = false;
    if(CMD.find("--chrsizes")== cmd_itr_end) { printf("chrsizes required. "); cmd_missing = true; }
    if(CMD.find("--snp")     == cmd_itr_end) { printf("snp      required. "); cmd_missing = true; } 
    if(CMD.find("--genome")  == cmd_itr_end) { printf("genome   required. "); cmd_missing = true; }
    if(CMD.find("--chrom")   == cmd_itr_end) { printf("chrom    required. "); cmd_missing = true; }
    if(CMD.find("--start")   == cmd_itr_end) { printf("start    required. "); cmd_missing = true; } 
    if(CMD.find("--end")     == cmd_itr_end) { printf("end      required. "); cmd_missing = true; } 
    if(CMD.find("--gff")     == cmd_itr_end) { printf("gff      required. "); cmd_missing = true; }
    if(CMD.find("--folder")  == cmd_itr_end) { printf("folder   required. "); cmd_missing = true; }
    if(vis_annotation==1 && fconsensus.size()==0) 
                                             { printf("setting of --consen with -vis required."); 
                                                                              cmd_missing = true; }
    if(cmd_missing)
    {
        printf("Exited.\n");
        exit(1);
    }
}
