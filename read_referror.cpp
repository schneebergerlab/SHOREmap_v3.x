// pls check the header file: read_referror.h for more info about this function. //
// this function is similar to read_chromosome(FILE* fp).
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include "globals.h"
#include "is_number.h"

bool read_referror(char* freferror)
{
    if(verbose) printf("Reading ref errors from file:\t\t%s...", freferror);
    
    FILE* fp  = fopen(freferror, "r");
    if(fp == NULL)
    {
        printf("Ref-error File \'%s\' does NOT exist. Exited.\n", freferror);
        exit(1);
    }
    if (!fp) return false;
    // caution size of word
    char word[1024];         
    std::string id_chr;
    std::string po_err; 
    // global map
    REFERROR.clear(); 
    while(!feof(fp))
    {
        id_chr = "";
        po_err = "";
        
        // read  chroms id
        fscanf(fp, "%s\n", word);
        // check chroms id
        std::string id_chr(word);
        if(CHR2SIZE.find(id_chr) == CHR2SIZE.end()) {//1.ignore entry not in given list
            fscanf(fp, "%s\n", word);
            continue;
        }
        
        // read  error pos
        fscanf(fp, "%s\n", word);  
        // check error pos
        if(!is_number(word))              continue; // 2.ignore entry that is not a number
        if(atol(word) > CHR2SIZE[id_chr]) continue; // 3.ignore entry that  exceeds bounds
        
        std::string po_err(word);
        // record if not found in map
        if(id_chr.length()>0 && po_err.length()>0)
            REFERROR.insert(std::pair<std::string,bool>(id_chr+".#."+po_err, true));
    }
    fclose(fp);
    
    // if no records found
    if (REFERROR.size()==0) 
    {
        if(CMD.find("-verbose")!=CMD.end()) 
            printf("no error reference provided: \n"); 
        return false;
    } 
    // print records
    if(verbose) printf("done.\n");
    if (verbose)
    {
        printf("ref errors provided. \n");
        map<std::string, bool>::iterator iter;
        for (iter = REFERROR.begin(); iter != REFERROR.end(); iter++) 
        {
          printf("\t\t\t\t\t%s, %ld\n", (*iter).first.c_str(), (unsigned long)(*iter).second);
        }
     }
    return true;
}
