// pls check the header file: print_allele_count.h for more info about this function. //
// CHR2POS2ALLELE1_COUNT, CHR2POS2ALLELE2_COUNT, and CHR2POS2ERROR_COUNT              //
/* Date  : 2013-March-08~09                                                           */     
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include "globals.h"

// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>

//#include "SHOREmap_plot.h"

bool print_allele_count()
{
    char out_file[1024];
    char pdf_file[1024];
    FILE* fp;
    
    DIR* dir_out_folder = opendir(out_folder.c_str());
    if(dir_out_folder  == NULL)
    {
        printf("Cannot create output path. Exited.\n");
        exit(1);
    }
    
    // part 1. write allele counts to file.
    sprintf(out_file, "%s/SHOREmap.winsize1.txt\0", out_folder.c_str());
    if(verbose) printf("Writing allele count info to file:\t%s...", out_file);
    fp = fopen(out_file, "w");
    if(!fp) {printf("ERROR: cannot open zoom info file: %s. Exited.\n", out_file); exit(1);}

    map<std::string, map<std::string, unsigned long> >::iterator itr4chr;
    for(itr4chr = CHR2POS2ALLELE1_COUNT.begin(); itr4chr != CHR2POS2ALLELE1_COUNT.end(); itr4chr ++)
    {
        std::string chrms = (*itr4chr).first;
        map<std::string, unsigned long> pos_cnt = CHR2POS2ALLELE1_COUNT[chrms];// position and count
        map<std::string, unsigned long>::iterator itr4pos;
        for(itr4pos=pos_cnt.begin(); itr4pos != pos_cnt.end(); itr4pos ++)
        {
            fprintf(fp, "%s\t", chrms.c_str());                         // chr id
            std::string postn = (*itr4pos).first;                       // position       
            fprintf(fp, "%s\t%ld\t", postn.c_str(), (*itr4pos).second); // count: allele 1, 2, error
            fprintf(fp, "%ld\t%ld\n",
                CHR2POS2ALLELE2_COUNT[chrms][postn], CHR2POS2ERROR_COUNT[chrms][postn]);
        }
    }
    if(!fp) fclose(fp);
    if(verbose) printf("done.\n");
    
    /* BELOW NOT USED - 2013-05-28 11:45 */
    // part 2. Carry out statistical analysis (calling R in the original implementation) 
    // sprintf(pdf_file, "%s/SHOREmap.pdf\0", out_folder.c_str());
    // if(verbose) printf("Plotting allele count info to file:\t%s...\n", pdf_file);
    /* call SHOREmap_plot.cpp: This funtion prepares data for plotting   */
    /* how to create and plot pdf file? this has to be resolved. Mar.09 - solved with Dislin  */
    
    closedir(dir_out_folder);
    
    return true;
}
