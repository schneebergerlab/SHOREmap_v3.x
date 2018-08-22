#include <fstream>
#include     <map>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>

#include "globals.h"

bool read_extract_cons(char* readfile)
{
    
    std::string line;
    unsigned long item        = 0;
    unsigned long con_mkr_num = 0; // number of con-bases related to markers
    
    std::ifstream FileCONSEN (readfile);
    if(!FileCONSEN.is_open())
    {
        printf("Consensus file \'%s\' does NOT exist. Exited.\n", readfile);
        exit(1);
    }
    if (verbose) printf("Reading base/error counts from file:\t%s...", readfile);
    
    /* prepare out file */
    char writefile[1024];  // caution
    if(out_folder[out_folder.length()-1] == '/')
    {
        sprintf(writefile, "%sextracted_consensus_%ld.txt\0", out_folder.c_str(), row_first);
    }
    else
    {
        sprintf(writefile, "%s/extracted_consensus_%ld.txt\0", out_folder.c_str(), row_first);
    }
    FILE* fp_write = fopen(writefile, "w");
    if(fp_write == NULL)  { printf("cannot open file to write consensus.\n"); exit(1);}

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
        pch  = strtok((char*)line.c_str(), token); // 1.chrid
        chr += pch;
        
        pch  = strtok(NULL, token);                // 2.pos
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
        /* // output selected columns
        if(firstline2) 
        {
            fprintf(fp_write, "%s\t%s\t", chr.c_str(), pos.c_str()); // chrid, position
            firstline2 = false;
        }
        else
        {
            fprintf(fp_write, "\n%s\t%s\t", chr.c_str(), pos.c_str()); // chrid, position
        }
        
        pch  = strtok(NULL, token);          // 3.base call
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 4. cov
        fprintf(fp_write, "%s\t", pch);
        
        pch   = strtok(NULL, token);         // 5. A_cnt
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 6.C_cnt
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 7.G_cnt
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 8.T.cnt
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 9.D.cnt
        fprintf(fp_write, "%s\t", pch);
        
        pch = strtok(NULL, token);           // 10.N.cnt
        fprintf(fp_write, "%s", pch);
                                             // 11. avg_hits
                                             // 12.avg_mm, 13., 14. ...
        */
    }
    if(verbose) 
    printf("\t%ld consen base info in total; %ld relate to markers provided.\n", item, con_mkr_num);
    
    if(!fp_write) fclose(fp_write);
    FileCONSEN.close();
    if(con_mkr_num <= 0)
    return false;
    else
    return true;
}
