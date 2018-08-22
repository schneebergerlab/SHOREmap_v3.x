// pls check the header file: write_log.h for more info about this function.                 //
// this function is required in the previous implementation, as it call R which requires the //
// 5 parameters including: reg_chromosome, reg_begin, reg_end, reg_freq_min, reg_freq_max.   //
// therefore, this function could be removed in this c-implementation.                       //
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include "globals.h"

// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>

bool write_zoom_region()
{
   char out_file[1024];
    
    DIR* dir_out_folder = opendir((char*)out_folder.c_str());
    if(dir_out_folder  == NULL)
    {
        printf("Cannot create output path. Exited.\n");
        exit(1);
    }   
    sprintf(out_file, "%s/SHOREmap_zoom_region.txt\0", (char*)out_folder.c_str());
    if(verbose) printf("Writing zoom info to file:\t\t%s...", out_file);
    
    FILE* fp = fopen(out_file, "w");
    if(!fp) {printf("ERROR: cannot open zoom info file: %s. Exited.\n", out_file); exit(1);}
    fprintf(fp, "%s \t %ld \t %ld \t %.2f \t %.2f \t\n", 
        reg_chromosome.c_str(), reg_begin, reg_end, reg_freq_min, reg_freq_max);

    fclose(fp);
    closedir(dir_out_folder);
    
    if(verbose) printf("done.\n");
    return true;
}
