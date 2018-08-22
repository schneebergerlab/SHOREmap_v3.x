#include     <map>
#include <cstring>
#include  <vector>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "globals.h"

bool read_allele_counts2(char* fconsensus)
{
    std::string line;
    unsigned long item = 0;
    
    std::ifstream FileCONSEN (fconsensus);
    if(!FileCONSEN.is_open())
    {
        printf("Consensus file \'%s\' does NOT exist. Exited.\n", fconsensus);
        exit(1);
    }
    if (verbose) printf("Reading consensus call info from file:\t%s: \n", fconsensus);

    bool firstline = true; 
    const char* token;
    while(FileCONSEN.good())
    {
        getline(FileCONSEN, line);
        if (line.length() == 0) {continue;}                                             // null line
        if(firstline)             // find out token of line. " " or "\t"; caution on for other cases
        {
            if(line.find(" ") != std::string::npos) { token = " ";}
            else if(line.find("\t") != std::string::npos) { token = "\t";}
            firstline = false;
        }
        item ++;
        if(item%10000000 == 0)  printf("    %ld consen base info read.\n", item);
        
        unsigned long cov   = 0;
        unsigned long A_sup = 0;
        unsigned long C_sup = 0;
        unsigned long G_sup = 0;
        unsigned long T_sup = 0;
        unsigned long D_sup = 0;
        unsigned long N_sup = 0;
        double      avg_hit = 0;
        std::string chr("");
        std::string pos("");
        std::string bas("");
        
        unsigned long ind = 0;
        char*  pch;
        pch  = strtok((char*)line.c_str(), token); // 1.chrid
        chr += pch;
        pch  = strtok(NULL, token);                // 2.pos
        pos += pch; 
        
        std::string ale_id = chr+".#."+pos;
        //if the position is not a recorded marker pos, ignore it (and related info).
        if(ALLELE1.find(ale_id) == ALLELE1.end()) {continue;}
        
        pch   = strtok(NULL, token);               // 3.base call
        bas  += pch;      
        pch   = strtok(NULL, token);               // 4. cov
        cov   = atoi(pch);   
        pch   = strtok(NULL, token);               // 5. A_cnt
        A_sup = atoi(pch);       
        pch   = strtok(NULL, token);               // 6.C_cnt
        C_sup = atoi(pch);     
        pch   = strtok(NULL, token);               // 7.G_cnt
        G_sup = atoi(pch);    
        pch   = strtok(NULL, token);               // 8.T.cnt
        T_sup = atoi(pch);  
        pch   = strtok(NULL, token);               // 9.D.cnt
        D_sup = atoi(pch); 
        pch   = strtok(NULL, token);               // 10.N.cnt
        N_sup = atoi(pch);
        pch   = strtok(NULL, token);               // 11.avg_hits
        avg_hit = atof(pch);
                                                   // 12.avg_mm, 13., 14. ...
                                                   
        // check fg_N, fg_LNDEL -- new on 2014-04-16 22:51
        if(D_sup>fg_INDEL_cov || N_sup>fg_N_cov || avg_hit>marker_hit) continue;                                           
        
        std::string allele1 = ALLELE1[ale_id];     // mut -- already swapped as original ref if background2 is 1 - af to calculate
        std::string allele2 = ALLELE2[ale_id];     // ref 
        
        unsigned long count_allele1 = 0;
        unsigned long count_allele2 = 0;
        unsigned long count_error   = 0;
        unsigned long coverage      = 0; // coverage = count_error + count_allele1 + count_allele2;
        unsigned long count_lower_allele  = 0;
        unsigned long count_higher_allele = 0;
        
        if (allele1 == "A")      count_allele1 = A_sup;
        else if (allele1 == "C") count_allele1 = C_sup;
        else if (allele1 == "G") count_allele1 = G_sup;
        else if (allele1 == "T") count_allele1 = T_sup;
        else if (allele1 == "-") count_allele1 = D_sup;
        else if (allele1 == "N") count_allele1 = N_sup;
        else ;
        
        if (allele2 == "A")      count_allele2 = A_sup;
        else if (allele2 == "C") count_allele2 = C_sup;
        else if (allele2 == "G") count_allele2 = G_sup;
        else if (allele2 == "T") count_allele2 = T_sup;
        else if (allele2 == "-") count_allele2 = D_sup; // not counted in original?
        else if (allele2 == "N") count_allele2 = N_sup; // not counted in original?
        else ;
        
        if(!(allele1 == "A") && !(allele2 == "A")) count_error += A_sup;
        if(!(allele1 == "C") && !(allele2 == "C")) count_error += C_sup;
        if(!(allele1 == "G") && !(allele2 == "G")) count_error += G_sup;
        if(!(allele1 == "T") && !(allele2 == "T")) count_error += T_sup;
        if(!(allele1 == "-") && !(allele2 == "-")) count_error += D_sup; // not counted in original?
        if(!(allele1 == "N") && !(allele2 == "N")) count_error += N_sup; // not counted in original?
        
        if(count_error < 0) {printf("%s\t%s error count < 0\n", chr.c_str(), pos.c_str()); exit(1);}
        
        coverage = count_error + count_allele1 + count_allele2; // caution on diff btw cov&coverage
        
        CHR2POS2ALLELE1_COUNT[chr][pos] = count_allele1;        // ref or pb allele
        CHR2POS2ALLELE2_COUNT[chr][pos] = count_allele2;        // mut or pa allele
        CHR2POS2ERROR_COUNT[chr][pos]   = count_error;          // err
        
        /* new 2013-06-13 23:37 */
        /* check if this SNP from parent B, we should swap the allele count for AF calc - no need, because alleles have been swapped when creating markers in this case   */
        /* QUALITY1 format is: FLAG4parentB#benji,quality */
        /*
        std::size_t found = QUALITY1[chr+".#."+pos].find("FLAG4parentB");
        if (found != std::string::npos)
        {
            // 2013-09-20 always interested in parent A - related allele frequency!!! - 2017-09-04
            CHR2POS2_ale1_ale2_err_COUNT[chr][atol(pos.c_str())] = (TRIPLE){count_allele1, count_allele2, count_error};
            
            // 2013-07-09 - no above-swap when recording consensus info, thus normal output in print_filtered_marker.cpp
            // CHR2POS2_ale1_ale2_err_COUNT[chr][atol(pos.c_str())] = (TRIPLE){count_allele1, count_allele2, count_error};
        }
        else
        {
            CHR2POS2_ale1_ale2_err_COUNT[chr][atol(pos.c_str())] = (TRIPLE){count_allele2, count_allele1, count_error};
        }
        */
        // 2017-09-04
        CHR2POS2_ale1_ale2_err_COUNT[chr][atol(pos.c_str())] = (TRIPLE){count_allele2, count_allele1, count_error};
        //
    }
    FileCONSEN.close();
    
    if(true) // output for checking: to turn off on 2014-04-13 19:36
    {
        std::string consensus_info_4_analysis = (out_folder + "SHOREmap_consensus_info_4_analysis.txt");
        fstream fout;
        fout.open (consensus_info_4_analysis.c_str(), ios::out);
        if(!fout.is_open())
        {
            printf("Cannot open file to write filtered markers (in print_filtered_marker.cpp). ");
            printf("Exited. \n");
            exit(1);
        }
        fout << "#chr\tpos\tallele1(af-to-calculate-if-background2=0)\tallele2\terror\tallele1\tallele2\tparent" << endl;
        
        map<std::string, map<unsigned long, TRIPLE> >::iterator it;
        for (it = CHR2POS2_ale1_ale2_err_COUNT.begin(); it!=CHR2POS2_ale1_ale2_err_COUNT.end(); it++)
        {
            std::string chrid = (*it).first;
            map<unsigned long, TRIPLE>::iterator it2;
            map<unsigned long, TRIPLE>::iterator it2_end;
            it2_end  = CHR2POS2_ale1_ale2_err_COUNT[chrid].end();
            for (it2 = CHR2POS2_ale1_ale2_err_COUNT[chrid].begin(); it2 != it2_end; it2 ++)
            {
                fout << chrid << "\t" << (*it2).first;
                
                std::stringstream ss;                    // position:    long       to        string
                ss << (*it2).first;
                std::string fndkey = (chrid+".#."+ss.str());
                /* // no need: because alleles have been swapped when creating markers.
                std::size_t found = QUALITY1[fndkey].find("FLAG4parentB");
                if (found != std::string::npos)
                {
                    fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[0]; // af to calculate:  
                    fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[1];
                }
                else
                {
                    fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[1]; // af to calculate:
                    fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[0];
                }
                */
                
                // 
                fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[1]; // af to calculate:
                fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[0];
                //
                fout << "\t" << CHR2POS2_ale1_ale2_err_COUNT[chrid][(*it2).first].Ci[2];
                
                fout << "\t" << ALLELE1[fndkey] << "\t"; // ALLELE1 - af to calculate
                fout << "\t" << ALLELE2[fndkey] << "\t";
                fout << "\t" << QUALITY1[fndkey] << endl;
            }
        }
        fout.close();
        printf("Consensus information has been recorded in SHOREmap_consensus_info_4_analysis.txt.\n");
    }
    if(verbose) 
    {
        printf("SUMMARY: %ld consensus base call info read; ", item);
        unsigned long num_marker_consen = 0;
        map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr;
        map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr_end;
        chr_itr     = CHR2POS2_ale1_ale2_err_COUNT.begin();
        chr_itr_end = CHR2POS2_ale1_ale2_err_COUNT.end();
        while(chr_itr != chr_itr_end)
        {
            num_marker_consen += (*chr_itr).second.size();
            chr_itr ++;
        }
        printf("%ld consensus info recorded for (part of) the provided markers. \n",
               num_marker_consen);
        if(INF>fg_INDEL_cov && INF>fg_N_cov && INF>marker_hit)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on indels, N and average hits.\n\n");
        }
        else 
        if(INF>fg_INDEL_cov && INF>fg_N_cov)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on indels, N.\n\n");
        }
        else
        if(INF>fg_INDEL_cov && INF>marker_hit)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on indels and average hits.\n\n");
        }
        else
        if(INF>fg_N_cov && INF>marker_hit)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on N and average hits.\n\n");
        }
        else
        if(INF>fg_INDEL_cov)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on indels.\n\n");
        }
        else
        if(INF>fg_N_cov)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on N.\n\n");
        }
        else
        if(INF>marker_hit)
        {
            printf("Consensus info (and thus markers to use) recorded must satisfy constaints on ");
            printf("coverage on average hits.\n\n");
        }
    }
    if(CHR2POS2ALLELE1_COUNT.size()<=0) 
    {
        return false;
    }
    else
    {
        return true;
    }
}
