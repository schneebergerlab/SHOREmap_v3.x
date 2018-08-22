/* NOTE: all snps related to each sample are obtained by aligning reads to a reference genome. 
   this function filters targeted SNPs of foreground (fg) according to background (bg):
      With two sequenced genomes, we use the targeted one as foreground, and the other as   
   background to reduce the number of SNPs for analysis. Specifically, suppose we have a        
   SNP M in fg A, we further check quality_reference (ref-base-call with high quality) 
   of bg B:                          
      If the reference allele of SNP M appears in quality_reference of bg B but with bad       
   support (too low or too high), quality,  or concordance, we might discard this SNP M (it cannot 
   be used as a marker in future analysis).            
      Otherwise, we keep SNP M as a marker for future analysis.                                    

2013-06-12 13:35 remove SNPs that do not have support in the bg-ref file, i.e., if cannot find
 the position of a SNP, remove this SNP???? -  have to verify effectiveness                       */

#include   <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include    <string>
#include   <sstream>
#include  <iostream>
#include    <vector>
#include       <map>

#include      "globals.h"
#include "split_string.h"
using namespace std;

int check_ref_err(char*  fg_snp_file,
                  char*  fg_consencall_file, 
                  char*  bg_ref_base_file, 
                  unsigned long min_cov,
                  unsigned long max_cov,
                  double        min_ccd, 
                  unsigned long max_N, 
                  unsigned long max_ID, 
                  unsigned long min_refq,
                  std::string* fg_snp_file_filtered,
                  map<std::string, MARK6> * fg_snp_map_filtered)
{
    // fg_snp_file          fg-SNP file
    // fg_consencall_file   fg-consensus info of SNPs in fg_snp_file (from SHOREmap extract)    
    // bg_ref_base_file     bg-ref file to help reduce the given fg-SNPs
    // min_cov              coverage                     of a bg-ref base
    // max_cov              coverage                     of a bg-ref base
    // min_ccd              concordance                  of a bg-ref base
    // min_refq             minimum quality              of a bg-ref base
    // max_N                maximum predicted ambiguity  of a fg-consensus base position
    // max_ID               maximum predicted '-'        of a fg-consensus base position    
    
    unsigned long SNP_in_bg_ref = 0;
    map<std::string, unsigned long> SNP_found_in_bg_ref;
    FILE* fp = fopen(fg_snp_file, "r");           //  e.g., SHOREmap_marker.bg_corrected_q25_f20_EMS
    if(!fp)
    {
        printf("Cannot open marker file. Exited.\n");
        exit(1);
    }
    unsigned long num_entry = 0;
    unsigned long num_raw   = 0;
    char pro[32];
    map<std::string, map<unsigned long, MARK6> > newMkrB4Map;            // < chr, <pos, {detail}> >
    while(!feof(fp))
    {
        // caution on file format: with exact 9 clolumns
        // char pro[32];
        char chr[32];
        char pos[32];
        char ref[32];
        char mut[32];
        char qua[32];
        char cov[32];
        char ccd[32];
        char aht[32];
        
        /* format required */
        if(!fscanf(fp, "%s\t", pro)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", chr)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", pos)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", ref)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", mut)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", qua)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", cov)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\t", ccd)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(fp, "%s\n", aht)) {printf("Error reading. exited.\n");exit(1);}
        
        /* added on 2013-06-14 09:54               */
        num_raw ++;
        if(atof(ccd) < reg_freq_min)        continue;// minimum frequency of fg-SNP
        if(atol(cov) < filter_min_coverage) continue;// minimum coverage  of fg-SNP
        if(atol(cov) > filter_max_coverage) continue;// maximum coverage  of fg-SNP
        if(atol(qua) < marker_score)        continue;// minimum score     of fg-SNP
        
        if(atol(qua) > quality_max) quality_max = atol(qua);
        
        /////////////////////
        std::string keychr = "";
        keychr += (std::string)chr;
        
        MARK6 other;
        other.ref  = "";
        other.ref += ref;
        other.mut  = "";
        other.mut += mut;
        other.qua  = "";
        other.qua += qua;
        other.cov  = "";
        other.cov += cov;
        other.ccd  = "";
        other.ccd += ccd;
        other.aht  = "";
        other.aht += aht;
        
        map<std::string, map<unsigned long, MARK6> >::iterator chr_itr;
        chr_itr = newMkrB4Map.find(keychr);
        if(chr_itr != newMkrB4Map.end())
        {
            (*chr_itr).second.insert(std::pair<unsigned long, MARK6>((unsigned long)atol(pos), other));
            num_entry ++;
        }
        else
        {
            map<unsigned long, MARK6> temp_mkr6;
            temp_mkr6.insert(std::pair<unsigned long, MARK6>((unsigned long)atol(pos), other));
            newMkrB4Map.insert(std::pair<std::string, map<unsigned long, MARK6> >(keychr, temp_mkr6));
            num_entry ++;
        }
    }
    fclose(fp);
    printf("%ld of (%ld) SNPs kept b4 checking background reference.\n", num_entry, num_raw);
    
    /* filter fg-SNPs with bg-ref bases. E.G., .../ConsensusAnalysis/quality_reference.txt        */
    FILE* reffp = fopen(bg_ref_base_file, "r");
    if(!reffp)
    {
        printf("Cannot open ref consensus file. exited.\n");
        exit(1);
    }
    while(!feof(reffp))
    {   // caution on file format: with exact 9 clolumns
        char pro2[32];
        char chr[32];
        char pos[32];
        char ref[32];
        char mut[32];
        char qua[32];
        char cov[32];
        char ccd[32];
        char aht[32];
        
        if(!fscanf(reffp, "%s\n", pro2)) {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", chr))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", pos))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", ref))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", mut))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", qua))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", cov))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", ccd))  {printf("Error reading. exited.\n");exit(1);}
        if(!fscanf(reffp, "%s\n", aht))  {printf("Error reading. exited.\n");exit(1);}
        
        if(atol(qua) > quality_max) quality_max = atol(qua);
        
        std::string chrkey = "";
        chrkey = (std::string)chr;
        
        map<std::string, map<unsigned long, MARK6> >::iterator chr_itr;
        chr_itr = newMkrB4Map.find(chrkey);
        if(chr_itr != newMkrB4Map.end())                                // same                  chr
        {
            map<unsigned long, MARK6>::iterator pos_itr;
            pos_itr = (*chr_itr).second.find((unsigned long)atol(pos));
            if(pos_itr != (*chr_itr).second.end())                      // same                  pos 
            {   
                if(strncmp(ref, (*pos_itr).second.ref.c_str(), 1) == 0) // same                  ref   
                {
                    bool is_erase = false;
                    if(atof(ccd)  < min_ccd)                            // min           concordance
                    {
                        is_erase  = true;
                    }
                    else if(atof(cov)<min_cov || atof(cov)>max_cov)     // min/max  coverage/support
                    {
                        is_erase  = true;
                    }
                    else if(atof(qua) < min_refq)                       // min      base     quality
                    {
                        is_erase  = true;
                    }                    
                    else 
                    {
                        ;
                    }
                    if(is_erase)
                    {
                        (*chr_itr).second.erase(pos_itr);
                        num_entry --;
                    }
                    else                             // record the bg-ref info for future clustering
                    {
                        std::string stemp = (string)chrkey+".#."+(string)pos;
                        std::string itemp = (string)qua+"#"+(string)cov+"#"+(string)ccd;
                        bgREF.insert(std::pair<string, string>(stemp, itemp)); // backcross
                    }
                    std::string fndkey = "";
                    fndkey += (*chr_itr).first;
                    fndkey += ".#.";
                    fndkey += (std::string)pos;
                    SNP_found_in_bg_ref.insert(std::pair<std::string, unsigned long>(fndkey, 1)); // all the supported positions
                }
            }
        }
    }
    fclose(reffp);
    printf("%ld SNPs after checking background reference file.\n", num_entry);
    
    /* remove a fg-snp without bg-ref support : 2013-06-12 15:02                                  */
    /* caution: this can result removal of a lot of targeted fg-SNPs.                             */
    /* TODO: not simply remove, instead, check cov, quality, ccd etc to decide to remove or not   */
    map<std::string, map<unsigned long, MARK6> >::iterator ichr_itr;
    map<std::string, map<unsigned long, MARK6> >::iterator ichr_itr_end;
    ichr_itr     = newMkrB4Map.begin();
    ichr_itr_end = newMkrB4Map.end();
    while(ichr_itr != ichr_itr_end)
    {
        map<unsigned long, MARK6>::iterator ipos_itr;
        map<unsigned long, MARK6>::iterator ipos_itr_end;
        ipos_itr     = (*ichr_itr).second.begin();
        ipos_itr_end = (*ichr_itr).second.end();
        while(ipos_itr != ipos_itr_end)
        {
            std::string tmpkey = "";
            tmpkey += (*ichr_itr).first;
            tmpkey += ".#.";
            stringstream ss;
            ss << (*ipos_itr).first;
            tmpkey += ss.str();
            // check: if a targeted fg-SNP does not have bachground ref-info, simply discard it.
            if(SNP_found_in_bg_ref.find(tmpkey) == SNP_found_in_bg_ref.end())
            {
                num_entry --;
                (*ichr_itr).second.erase(ipos_itr++);
            }
            else            
            ipos_itr ++;
        }
        ichr_itr ++;
    }
    printf("%ld SNPs after removing those not supported with bg-reference.\n", num_entry);
    
    /* the following filtering is carried out only when fg-consensus information is provided      */
    /* check foreground consensus call info file. E.G., .../myVisual/extracted_consensus_0408.txt */
    /* using number of '-' and/or 'N' if their default values are re-set by user                  */
    if((strlen(fg_consencall_file) != 0) && (max_N != INF || max_ID != INF))
    {
        FILE* ownCfp = fopen(fg_consencall_file, "r");
        if(!ownCfp)
        {
            printf("Cannot open ref consensus file for checking coverage of 'N'/'-'. exited.\n");
            exit(1);
        }
        char word[1024];     // caution
        while(!feof(ownCfp))
        {
            // caution on file format: with 65 clolumns, while we are interested in the first 10 columns 
            unsigned long cov   = 0;
            unsigned long A_sup = 0;
            unsigned long C_sup = 0;
            unsigned long G_sup = 0;
            unsigned long T_sup = 0;
            unsigned long D_sup = 0;
            unsigned long N_sup = 0;
        
            std::string chr("");
            std::string pos("");
            std::string bas("");
        
            // chromosome id
            if(!fscanf(ownCfp, "%s\n", word)) {printf("Error reading. exited.\n");exit(1);}
            chr += word;  
            // base position
            if(!fscanf(ownCfp, "%s\n", word)) {printf("Error reading. exited.\n");exit(1);}
            pos += word;
            
            std::string chrkey = "";
            chrkey = (std::string)chr;
            
            map<std::string, map<unsigned long, MARK6> >::iterator chr_itr;
            chr_itr = newMkrB4Map.find(chrkey);
            if(chr_itr != newMkrB4Map.end())                                // same chr
            {
                map<unsigned long, MARK6>::iterator pos_itr;
                pos_itr = (*chr_itr).second.find((unsigned long)atol(pos.c_str()));
                if(pos_itr != (*chr_itr).second.end())                      // same pos 
                {
                    // note: this coverage = 'A'+'C'+'G'+'T'+'-'+'N', '-': deletion, 'N': ambiguition
                    if(!fscanf(ownCfp, "%s\n", word)) {printf("Error reading. exited.\n");exit(1);}
                    cov= atoi(word); 
            
                    if(!fscanf(ownCfp, "%ld\n", &A_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(!fscanf(ownCfp, "%ld\n", &C_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(!fscanf(ownCfp, "%ld\n", &G_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(!fscanf(ownCfp, "%ld\n", &T_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(!fscanf(ownCfp, "%ld\n", &D_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(!fscanf(ownCfp, "%ld\n", &N_sup)) {printf("Error reading. exited.\n");exit(1);}
                    if(N_sup>max_N || D_sup>max_ID)
                    {
                        (*chr_itr).second.erase(pos_itr);
                        num_entry --;
                    }
                }
            }
            // ignore remaining info
            while(getc(ownCfp)!='\n' && !feof(ownCfp)); // caution
        }
        fclose(ownCfp);
        printf("%ld SNPs after checking '-'/'N' consensus file.\n", num_entry);
    }
    else
    {
        printf("SNPs are not filtered according to N-call or indel-call on their positions.\n");
    }
    /* write retained targeted SNPs from fg_snp_file to a new file                                */
    char wfile[1024];
    // background reference base info are recorded in file name
    if(max_N<INF && max_ID<INF)
    {
        sprintf(wfile, "%sSHOREmap_filtered_on_bgREF_f%.3f_q%d_minc%ld_maxc%ld_maxN%ld_maxID%ld.txt\0", 
            out_folder.c_str(), min_ccd, (int)min_refq, min_cov, max_cov, max_N, max_ID);        
    }
    else if(max_N<INF)
    { 
        sprintf(wfile, "%sSHOREmap_filtered_on_bgREF_f%.3f_q%d_minc%ld_maxc%ld_maxN%ld.txt\0", 
            out_folder.c_str(), min_ccd, (int)min_refq, min_cov, max_cov, max_N);            
    }
    else if(max_ID<INF)
    {
        sprintf(wfile, "%sSHOREmap_filtered_on_bgREF_f%.3f_q%d_minc%ld_maxc%ld_maxID%ld.txt\0", 
            out_folder.c_str(), min_ccd, (int)min_refq, min_cov, max_cov, max_ID);  
    }
    else
    {
        sprintf(wfile, "%sSHOREmap_filtered_on_bgREF_f%.3f_q%d_minc%ld_maxc%ld.txt\0", 
            out_folder.c_str(), min_ccd, (int)min_refq, min_cov, max_cov);
    }
    *fg_snp_file_filtered += wfile;
    FILE* wfp = fopen(wfile, "w");
    map<std::string, map<unsigned long, MARK6> >::iterator chr_itr;
    chr_itr = newMkrB4Map.begin();
    while(chr_itr != newMkrB4Map.end())
    {
        map<unsigned long, MARK6>::iterator pos_itr;
        pos_itr = (*chr_itr).second.begin();
        while(pos_itr != (*chr_itr).second.end())
        {
            /* record the SNPs in file for the purpose of future check                            */
            fprintf(wfp, "%s\t%s\t%ld\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    pro, 
                    (*chr_itr).first.c_str(),
                    (*pos_itr).first,
                    (*pos_itr).second.ref.c_str(),
                    (*pos_itr).second.mut.c_str(),
                    (*pos_itr).second.qua.c_str(),
                    (*pos_itr).second.cov.c_str(),
                    (*pos_itr).second.ccd.c_str(),
                    (*pos_itr).second.aht.c_str());
            /* record the SNPs in returned map                                                    */
            /* search key */
            std::string chrpos = "";
            chrpos += (*chr_itr).first;
            chrpos += ".#.";
            std::stringstream ss;  // position: long to string
            ss << (*pos_itr).first;
            chrpos += ss.str();
            (*fg_snp_map_filtered).insert(std::pair<std::string, MARK6>(chrpos, (*pos_itr).second));
            /* next SNP */
            pos_itr ++; 
        }
        chr_itr ++;
    }
    fclose(wfp);

    return 1;
}
/* 
format of "filtered_snps_bgREF_concordance%.3f_quality%d.txt\0":
1.project-name 2.chrID 3.position 4.ref-base 5.mut-base 6.quality 7.coverage 8.concordance 9.avg_hit
e.g.:
Pep2-1         Aa.chr1 38610      T          T          40        119        1             1
*/
