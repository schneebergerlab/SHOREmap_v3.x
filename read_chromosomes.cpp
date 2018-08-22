// pls check the header file: read_chromosomes.h for more info about this function.               //
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include "globals.h"

#include "print_error_exit.h"

bool read_chromosomes(char* fchrsizes)
{
    FILE* fp = fopen(fchrsizes, "r");
    if (!fp) return false;  
    /* global max value of chr sizes */
    chrsizes_max = 0; 
    char word[1024];                                                         // caution size of word
    std::string id_chr;
    unsigned long value_chr = 0;
    /* global map */
    CHR2SIZE.clear();
    while(!feof(fp))
    {
        /*   id */
        id_chr = "";
        if(!fscanf(fp, "%s\n", word)) print_error_exit((char*)"read_chromosomes(...)", true);      
        std::string id_chr(word);
        /* value */
        if(!fscanf(fp, "%s\n", word)) print_error_exit((char*)"read_chromosomes(...)", true);     
        value_chr = atol(word);
        /* record if not found in map */
        if(CHR2SIZE.find(id_chr) == CHR2SIZE.end() && value_chr>0) 
        {
            CHR2SIZE.insert(std::pair<string,unsigned long>(id_chr, value_chr));
            if(value_chr > chrsizes_max) chrsizes_max = value_chr;
        }
        else 
        {
            printf("repeated chromosome id or invalid chromosome size detected.\n");
        }
    }
    fclose(fp);
    if (verbose)
    {
        if(CHR2SIZE.size() <= 10)
        {
            map<std::string, unsigned long>::iterator iter;
            for (iter = CHR2SIZE.begin(); iter != CHR2SIZE.end(); iter++) 
            {
                printf("\t\t\t\t\t%s, %ld\n", (*iter).first.c_str(), (*iter).second);
            }
        }
        else
        {
            printf("\t\t\t\t\tSize info recorded for %ld chromosomes. \n", CHR2SIZE.size());
        }
    }
    if (CHR2SIZE.size()==0) {return false;}
    return true;
}
