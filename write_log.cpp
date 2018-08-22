// pls check the header file: write_log.h for more info about this function. //
#include  <stdio.h>
#include   <string>
#include <stdlib.h>
#include "globals.h"

// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>

bool write_log()
{
    /* TODO: make it general to other sub-main functions  */
    char out_file[1024];
    
    DIR* dir_out_folder = opendir((char*)out_folder.c_str());
    if(dir_out_folder  == NULL)
    {
        printf("Cannot create output path. Exited.\n");
        exit(1);
    }
    
    sprintf(out_file, "%s/SHOREmap.log\0", (char*)out_folder.c_str());
    
    if(verbose) printf("Writing CMD log to file:\t\t%s...", out_file);
    
    FILE* fp = fopen(out_file, "w");
    if(!fp) {printf("ERROR: cannot open output file: %s. Exited.\n", out_file); exit(1);}
    
    fprintf(fp, "Mandatory: \n");
    if(CMD.find("--chrsizes")!=CMD.end())            
        fprintf(fp, " --chrsizes\t\t%s\n", CMD["--chrsizes"].c_str());                                  
    if(CMD.find("--folder")!=CMD.end())              
        fprintf(fp, " --folder\t\t%s\n", CMD["--folder"].c_str());                                  
    if(CMD.find("--marker")!=CMD.end())              
        fprintf(fp, " --marker\t\t%s\n", CMD["--marker"].c_str());                                 
    if(CMD.find("--consen")!=CMD.end())              
        fprintf(fp, " --consen\t\t%s\n", CMD["--consen"].c_str());
    
    fprintf(fp, "Other: \n");
    if      (CMD.find("-target")!=CMD.end())               
        fprintf(fp, "  -target\t\t%s\n", CMD["-target"].c_str());
    else                                             
        fprintf(fp, "  -target\t\t1.0 [default]\n");
    if      (CMD.find("-conf-int")!=CMD.end())             
        fprintf(fp, "  -conf-int\t\t%s\n", CMD["-conf-int"].c_str());
    else                                             
        fprintf(fp, "  -conf-int\t\t2 [default]\n");
    if     (CMD.find("--conf")!=CMD.end())                
        fprintf(fp, " --conf\t\t\t%s\n", CMD["--conf"].c_str());
    else                                             
        fprintf(fp, " --conf\t\t\t  [default]\n");
                                
    if     (CMD.find("--window-size")!=CMD.end())         
        fprintf(fp, " --window-size\t\t%s\n", CMD["--window-size"].c_str());
    else                                             
        fprintf(fp, " --window-size\t\t50000 [default]\n");
    if     (CMD.find("--window-step")!=CMD.end())        
        fprintf(fp, " --window-step\t\t%s\n", CMD["--window-step"].c_str());
    else                                             
        fprintf(fp, " --window-step\t\t10000 [default]\n");
    if     (CMD.find("--peak-window-size")!=CMD.end())    
        fprintf(fp, " --peak-window-size\t%s\n", CMD["--peak-window-size"].c_str());
    else                                         
        fprintf(fp, " --peak-window-size\t50000 [default]\n");
    if     (CMD.find("--peak-window-step")!=CMD.end())
        fprintf(fp, " --peak-window-step\t%s\n", CMD["--peak-window-step"].c_str());
    else                                         
        fprintf(fp, " --peak-window-step\t10000 [default]\n");

    if     (CMD.find("--min-marker")!=CMD.end())   	 
        fprintf(fp, " --min-marker\t\t%s\n", CMD["--min-marker"].c_str());
    else                                         
        fprintf(fp, " --min-marker\t\t0 [default]\n");
    if     (CMD.find("--min-coverage")!=CMD.end())    
        fprintf(fp, " --min-coverage\t\t%s\n", CMD["--min-coverage"].c_str());
    else                                         
        fprintf(fp, " --min-coverage\t\t0 [default]\n");
    if     (CMD.find("--max-coverage")!=CMD.end())    
        fprintf(fp, " --max-coverage\t\t%s\n", CMD["--max-coverage"].c_str());
    else                                         
        fprintf(fp, " --max-coverage\t\tINF [default]\n");

    if     (CMD.find("--outlier-window-size")!=CMD.end())
        fprintf(fp, " --outlier-window-size\t%s\n", CMD["--outlier-window-size"].c_str());
    else                                         
        fprintf(fp, " --outlier-window-size\t200000 [default]\n");
    if      (CMD.find("--outlier-pvalue")!=CMD.end())  
        fprintf(fp, " --outlier-pvalue\t%s\n", CMD["--outlier-pvalue"].c_str());
    else                                         
        fprintf(fp, " --outlier-pvalue\t0.05 [default]\n");
    if     (CMD.find("--mis-phenotyped")!=CMD.end())  
        fprintf(fp, " --mis-phenotyped\t%s\n", CMD["--mis-phenotyped"].c_str());
    else                                         
        fprintf(fp, " --mis-phenotyped\t0.00 [default]\n");

    if     (CMD.find("--chromosome")!=CMD.end())      
        fprintf(fp, " --chromosome\t\t%s\n", CMD["--chromosome"].c_str());
    else                                         
        fprintf(fp, " --chromosome\t\tNULL [default]\n");
    if     (CMD.find("--begin")!=CMD.end())           
        fprintf(fp, " --begin\t\t%s\n", CMD["--begin"].c_str());
    else                                         
        fprintf(fp, " --begin\t\t-1 [default]\n");
    if     (CMD.find("--end")!=CMD.end())             
        fprintf(fp, " --end\t\t\t%s\n", CMD["--end"].c_str());
    else                                         
        fprintf(fp, " --end\t\t\t-1 [default]\n");
    if      (CMD.find("--minfreq")!=CMD.end())         
        fprintf(fp, " --minfreq\t\t%s\n", CMD["--minfreq"].c_str());
    else                                         
        fprintf(fp, " --minfreq\t\t-1 [default]\n");
    if     (CMD.find("--maxfreq")!=CMD.end())         
        fprintf(fp, " --maxfreq\t\t%s\n", CMD["--maxfreq"].c_str());
    else                                         
        fprintf(fp, " --maxfreq\t\t-1 [default]\n");

    if     (CMD.find("--referrors")!=CMD.end())       
        fprintf(fp, " --referrors\t\t%s\n", CMD["--referrors"].c_str());
    else                                          
        fprintf(fp, " --referrors\t\tNULL [default]\n");
    if      (CMD.find("-background2")!=CMD.end())      
        fprintf(fp, "  -background2\t\t%s\n", CMD["-background2"].c_str());
    else                                         
        fprintf(fp, "  -background2\t\t0 [default]\n");
    if      (CMD.find("-verbose")!=CMD.end())          
        fprintf(fp, "  -verbose\t\t%s\n", CMD["-verbose"].c_str());
    else                                         
        fprintf(fp, "  -verbose\t\t0 [default]\n");
    
    if      (CMD.find("-boost-max")!=CMD.end())        
        fprintf(fp, "  -boost-max\t\t%s\n", CMD["-boost-max"].c_str());
    else                                         
        fprintf(fp, "  -boost-max\t\t10000 [default]\n");
    if      (CMD.find("-r-max")!=CMD.end())            
        fprintf(fp, "  -r-max\t\t%s\n", CMD["-r-max"].c_str());
    else                                         
        fprintf(fp, "  -r-max\t\t10000 [default]\n");
    // related to "-plot-r"
    if      (CMD.find("-plot-r")!=CMD.end())           
        fprintf(fp, "  -plot-r\t\t%s\n", CMD["-plot-r"].c_str());
    else                                         
        fprintf(fp, "  -plot-r\t\t0 [default]\n");
    if      (CMD.find("-plot-marker")!=CMD.end())      
        fprintf(fp, "  -plot-marker\t\t%s\n", CMD["-plot-marker"].c_str());
    else                                         
        fprintf(fp, "  -plot-marker\t\t0 [default]\n");
                             
    if      (CMD.find("-runid")!=CMD.end())            
        fprintf(fp, "  -runid\t\t%s\n", CMD["-runid"].c_str());
    else                                         
        fprintf(fp, "  -runid\t\t1 [default]\n");
    
    fclose(fp);
    closedir(dir_out_folder);
    
    if(verbose) printf("done.\n");
    return true;
}
