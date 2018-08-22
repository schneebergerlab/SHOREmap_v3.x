/* filter foreground markers according read support, coverage and allele concordance */
/* Date: 2013-03-26   */

#include      <stdlib.h>
#include       <stdio.h>
#include        <math.h>
#include       <ctype.h>
#include        <time.h>
#include           <map>
#include       <cstring>
#include       <ctype.h>

// file, dir
//#include    <dirent.h>
//#include  <sys/stat.h>
//#include <sys/types.h>
//#include    <unistd.h>

#include "globals.h"

int filter_fg(char* rfile, std::string* fout1, std::string* fout2, std::string* fout3)
{
    /* 1. read markers (possibly background-corrected) */
    FILE* cFileFG = fopen(rfile, "r");
    if(!cFileFG)
    {
        printf("Foreground file \'%s\' does NOT exist. Exited.\n", rfile);
        exit(1);
    }
    /* 2. prepare output files */
    char out_file1[1024];
    char out_file2[1024];
    char EMS_file3[1024];
    if(CMD.find("--bg") != CMD.end())
    {
       sprintf(out_file1, "%s", rfile);
    }
    else
    {
       sprintf(out_file1, "%s%s", out_folder.c_str(), "SHOREmap_marker.no_correction");
    }
    
    // new on 2013-10-05 14:35
    if(marker_hit < (double)INF)
    {
       sprintf(out_file1, "%s_mh%.4f\0", out_file1, marker_hit);
    }
    if(filter_min_coverage>0 && filter_max_coverage<INF)
    {
       sprintf(out_file1, "%s_ic%ld_ac%ld_q%ld\0", out_file1, filter_min_coverage, filter_max_coverage, marker_score);
       sprintf(out_file2, "%s_f%2.1f\0",           out_file1, marker_freq*100);
       sprintf(EMS_file3, "%s_EMS\0",              out_file2);
    }
    else if(filter_min_coverage>0)
    {
       sprintf(out_file1, "%s_ic%ld_q%ld\0",       out_file1, filter_min_coverage, marker_score);
       sprintf(out_file2, "%s_f%2.1f\0",           out_file1, marker_freq*100);
       sprintf(EMS_file3, "%s_EMS\0",              out_file2);
    }
    else if(filter_max_coverage<INF)
    {
       sprintf(out_file1, "%s_ac%ld_q%ld\0",       out_file1, filter_max_coverage, marker_score);
       sprintf(out_file2, "%s_f%2.1f\0",           out_file1, marker_freq*100);
       sprintf(EMS_file3, "%s_EMS\0",              out_file2);
    }
    else
    {
       sprintf(out_file1, "%s_q%ld\0",             out_file1, marker_score);
       sprintf(out_file2, "%s_f%2.1f\0",           out_file1, marker_freq*100);    // caution
       sprintf(EMS_file3, "%s_EMS\0",              out_file2);
    }
    /* return name of output files */
    *fout1 += (std::string)out_file1;
    *fout2 += (std::string)out_file2;
    *fout3 += (std::string)EMS_file3;
    FILE* fp_of1 = fopen(out_file1, "w");
    FILE* fp_of2 = fopen(out_file2, "w");
    if(!fp_of1 || !fp_of2)
    {
        printf("Cannot open file to write:\t%s, \n\t%s.\nExited.\n", out_file1, out_file2);
        exit(1);
    }
    if(verbose) 
    {   
        printf("\nFiltered SNPs will be written to the following files: \n");
        printf("\t%s\n", out_file1);
        printf("\t%s\n", out_file2);
        if(only_EMS == 1)
        printf("\t%s\n", EMS_file3);
    }
    /* 3. filter markers according to marker score and concordance */
    bool firstline1 = true;
    bool firstline2 = true;
    while(!feof(cFileFG))
    {          
        char projname[32];
        char chrid[32];
        char position[32]; // caution, not long.
        char refb;
        char mutb;
        char qscore[32];
        char pcoverage[32];
        char pconcordance[32];
        char pavghits[32];
        
        fscanf(cFileFG, "%s\n", projname);
        fscanf(cFileFG, "%s\n", chrid);
        fscanf(cFileFG, "%s\n", position);
        fscanf(cFileFG, "%c\n", &refb);
        fscanf(cFileFG, "%c\n", &mutb);
        fscanf(cFileFG, "%s\n", qscore);
        fscanf(cFileFG, "%s\n", pcoverage);
        fscanf(cFileFG, "%s\n", pconcordance);
        fscanf(cFileFG, "%s\n", pavghits);
        
        if(atof(pavghits) > marker_hit) continue; // caution: filtering with avg_hit
        if(filter_min_coverage>0 && atof(pcoverage)<filter_min_coverage)   continue; // filtering1.0
        if(filter_max_coverage<INF && atof(pcoverage)>filter_max_coverage) continue; // filtering1.1
        
        if(atol(qscore) >= marker_score)                                      // filtering2     
        {
            if(firstline1)
            {
                firstline1 = false;
            }
            else
            {
                fprintf(fp_of1, "\n");
            }
            fprintf(fp_of1, "%s\t%s\t%s\t%c\t%c\t%s\t%s\t%s\t%s", 
                projname, chrid, position, refb, mutb, qscore, pcoverage, pconcordance, pavghits);
        }
        if(atol(qscore) >= marker_score && atof(pconcordance) >= marker_freq) // filtering3
        {
            if(firstline2)
            {
                firstline2 = false;
            }
            else
            {
                fprintf(fp_of2, "\n");
            }
            fprintf(fp_of2, "%s\t%s\t%s\t%c\t%c\t%s\t%s\t%s\t%s", 
                projname, chrid, position, refb, mutb, qscore, pcoverage, pconcordance, pavghits);
        }
    }   
    fclose(fp_of1);
    fclose(fp_of2);
    fclose(cFileFG);
    
    /*4. retain only EMS-induced markers in out_file2 */
    /*   TODO: combined into the above process        */
    if(only_EMS==1)                                                         // filtering4
    {
        FILE* fp_of3 = fopen(EMS_file3, "w");
        if(!fp_of3)
        {
            printf("cannot open file to write info: \t%s\n", EMS_file3); exit(1);
        }
        FILE* OUT2 = fopen(out_file2, "r");
        if(!OUT2)
        {
            printf("File \'%s\' does NOT exist. Exited.\n", out_file2); exit(1);
        }
        firstline1 = true;
        long ln = 0;
        while(!feof(OUT2))
        {              
            char projname[32];
            char chrid[32];
            char position[32]; // caution: not long.
            char refb;
            char mutb;
            char qscore[32];
            char pcoverage[32];
            char pconcordance[32];
            char pavghits[32];
              
            fscanf(OUT2, "%s\n", projname);
            fscanf(OUT2, "%s\n", chrid);
            fscanf(OUT2, "%s\n", position);
            fscanf(OUT2, "%c\n", &refb);
            fscanf(OUT2, "%c\n", &mutb);
            fscanf(OUT2, "%s\n", qscore);
            fscanf(OUT2, "%s\n", pcoverage);
            fscanf(OUT2, "%s\n", pconcordance);
            fscanf(OUT2, "%s\n", pavghits);
            
            if((refb=='G' && mutb=='A') || (refb=='C' && mutb=='T')) // other "EMS-like" interested?
            {
                if(firstline1)
                {
                    firstline1 = false;
                }
                else
                {
                    fprintf(fp_of3, "\n");
                }
                fprintf(fp_of3, "%s\t%s\t%s\t%c\t%c\t%s\t%s\t%s\t%s",
                                projname, chrid, position, refb, mutb,
                                qscore, pcoverage, pconcordance, pavghits);
            }
        }
        fclose(OUT2);
        fclose(fp_of3);
    }  
    return 1;
}
