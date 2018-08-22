// pls check the header file: read_allele_counts.h for more info about this function. //
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include <iostream>
#include "globals.h"

bool read_allele_counts(char* fconsensus)
{
    /* TODO: replace this with read_allele_counts2(...) */
    FILE* fp = fopen(fconsensus, "r");
    if(fp == NULL)
    {
        printf("Consensus file \'%s\' does NOT exist. Exited.\n", fconsensus);
        exit(1);
    }
    if (verbose) printf("Reading base/error counts from file:\t%s...", fconsensus);
    
    char word[1024];     // caution
    while(!feof(fp))
    {
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
        fscanf(fp, "%s\n", word);
        chr += word;  
        // base position
        fscanf(fp, "%s\n", word);
        pos += word;
        
        std::string ale_id = chr+".#."+pos;
        // if the position is not a recorded marker pos, ignore it (and related info).
        if(ALLELE1.find(ale_id) == ALLELE1.end()) {while(getc(fp)!='\n' && !feof(fp)); continue;}
        
        //skip bas
        fscanf(fp, "%s\n", word);
        
        // coverage
        fscanf(fp, "%s\n", word);
        cov= atoi(word); // note: this cov = 'A'+'C'+'G'+'T'+'-'+'N', '-': deletion, 'N': ambiguious
        //unsigned long coverage = cov;
        
        // support of 'A'
        fscanf(fp, "%ld\n", &A_sup);
        // support of 'C'
        fscanf(fp, "%ld\n", &C_sup);
        // support of 'G'
        fscanf(fp, "%ld\n", &G_sup);
        // support of 'T'
        fscanf(fp, "%ld\n", &T_sup);
        // support of '-': deletion
        fscanf(fp, "%ld\n", &D_sup);
        // support of N:   no base call
        fscanf(fp, "%ld\n", &N_sup);
        
        std::string allele1 = ALLELE1[ale_id];
        std::string allele2 = ALLELE2[ale_id];
        
        unsigned long count_allele1 = 0;
        unsigned long count_allele2 = 0;
        unsigned long count_error   = 0;
        unsigned long coverage      = 0; // note: coverage = 'A'+'C'+'G'+'T' = 
                                         // count_error + count_allele1 + count_allele2;
        
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
        else if (allele2 == "-") count_allele2 = D_sup; // not counted in original script?
        else if (allele2 == "N") count_allele2 = N_sup; // not counted in original script?
        else ;
        
        if(!(allele1 == "A") && !(allele2 == "A")) count_error += A_sup;
        if(!(allele1 == "C") && !(allele2 == "C")) count_error += C_sup;
        if(!(allele1 == "G") && !(allele2 == "G")) count_error += G_sup;
        if(!(allele1 == "T") && !(allele2 == "T")) count_error += T_sup;
        if(!(allele1 == "-") && !(allele2 == "-")) count_error += D_sup; // not counted in original?
        if(!(allele1 == "N") && !(allele2 == "N")) count_error += N_sup; // not counted in original?
        
        coverage = count_error + count_allele1 + count_allele2; // caution on diff btwn cov&coverage
        
        CHR2POS2ALLELE1_COUNT[chr][pos] = count_allele1;
        CHR2POS2ALLELE2_COUNT[chr][pos] = count_allele2;
        CHR2POS2ERROR_COUNT[chr][pos]   = count_error;      
        
        // shq-added on 2013-03-10: [string_chr][long_pos] <--> TRIPLE_{allele1, allele2, error} //
        // this variable will play as 'data' in SHOREmap_plot of previous implementation         //
        CHR2POS2_ale1_ale2_err_COUNT[chr][atol(pos.c_str())] = (TRIPLE){count_allele1, count_allele2, count_error};
        
        // ignore remaining info
        while(getc(fp)!='\n' && !feof(fp));
    }
     
    fclose(fp);
    
    if (verbose) printf("done.\n");
    // check resutls
    if(verbose)
    {
        cout << "chr.#.pos\tallele1\tallele2\terror\tallele1\tallele2";
        cout << " (allele1 is usually the ref base, allele2 is the mut base)" << endl;
        map<std::string, map<std::string, unsigned long> >::iterator it;
        for (it = CHR2POS2ALLELE1_COUNT.begin(); it!=CHR2POS2ALLELE1_COUNT.end(); it++)
        {
            std::string chrid = (*it).first;
            map<std::string, unsigned long>::iterator it2;
            map<std::string, unsigned long>::iterator it2_end;
            it2_end = CHR2POS2ALLELE1_COUNT[chrid].end();
            for (it2 = CHR2POS2ALLELE1_COUNT[chrid].begin(); it2 != it2_end; it2++)
            {
                cout << chrid << ".#." << (*it2).first << "\t" << (*it2).second;
                cout << "\t" << CHR2POS2ALLELE2_COUNT[chrid][(*it2).first];
                cout << "\t" << CHR2POS2ERROR_COUNT[chrid][(*it2).first];
                cout << "\t" << ALLELE1[chrid+".#."+(*it2).first];
                cout << "\t" << ALLELE2[chrid+".#."+(*it2).first] <<endl;
            }
        }
    }
    if (CHR2POS2ALLELE1_COUNT.size()<=0) return false;
    return true;
}
