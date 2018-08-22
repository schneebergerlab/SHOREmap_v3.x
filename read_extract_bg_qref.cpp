#include <fstream>
#include     <map>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>

#include "globals.h"

/* read and extract quality_reference of background mutant for markers of foreground mutant       */
/* the extracted info will be used for filtering markers of foreground mutant later               */
bool read_extract_bg_qref(char* readfile)
{
    std::string line;
    unsigned long item        = 0;
    unsigned long con_mkr_num = 0; // number of con-bases related to markers
    
    std::ifstream FileCONSEN (readfile);
    if(!FileCONSEN.is_open())
    {
        printf("Background ref-base file \'%s\' does NOT exist. Exited.\n", readfile);
        exit(1);
    }
    if (verbose) printf("Reading ref-base info from file:\t%s...", readfile);
    
    /* prepare out file */
    char writefile[1024];  // caution
    if(out_folder[out_folder.length()-1] == '/')
    {
        sprintf(writefile, "%sextracted_quality_ref_base_%ld.txt\0", out_folder.c_str(), row_first);
    }
    else
    {
        sprintf(writefile, "%s/extracted_quality_ref_base_%ld.txt\0", out_folder.c_str(), row_first);
    }
    printf("Extracted quality reference in \t %s.\n", writefile);
    FILE* fp_write = fopen(writefile, "w");
    if(fp_write == NULL)  { printf("cannot open file to write ref-base info.\n"); exit(1);}

    bool firstline  = true; 
    bool firstline2 = true;
    const char* token;
    while(FileCONSEN.good())
    {
        getline(FileCONSEN, line);
        if (line.length() == 0) {continue;} // null line
        if(firstline)// separating token of line. " " or "\t"; for other cases, tune it specifically.
        {
            if(line.find(" ") != std::string::npos) { token = " ";}
            else if(line.find("\t") != std::string::npos) { token = "\t";}
            firstline = false;
        }
        item ++;
        std::string templine("");
        templine += line;
        
        std::string chr("");
        std::string pos("");
        std::string bas("");
        
        unsigned long ind = 0;
        char* pch;
        pch  = strtok((char*)line.c_str(), token); // 1.project  name
        
        pch  = strtok(NULL, token);                // 2.chromosome id
        chr += pch; 
        
        pch  = strtok(NULL, token);                // 3.snp  position
        pos += pch; 
        
        std::string ale_id = chr+".#."+pos;
        //if the position is not a recorded marker pos, ignore it (and related info).
        if(ALLELE1.find(ale_id) == ALLELE1.end()) {continue;}
        
        con_mkr_num ++;
        // output the whole line
        if(firstline2)
        {
            fprintf(fp_write, "%s\n", templine.c_str()); // chrid, position
            firstline2 = false;
        }
        else
        {
            fprintf(fp_write, "%s\n", templine.c_str()); // chrid, position
        }
    }
    if(verbose) 
    printf("    %ld bg-ref base info in total; %ld relate to markers provided.\n", item, con_mkr_num);
    
    if(!fp_write) fclose(fp_write);
    FileCONSEN.close();
    if(con_mkr_num <= 0) return false;
    else
    return true;
}
