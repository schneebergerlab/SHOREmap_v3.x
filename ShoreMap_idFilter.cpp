/* this function filters indels from two samples; output specific indels as well as common ones */
#include  <stdio.h>
#include <string.h>
#include <stdlib.h>
#include  <fstream>
#include <iostream>
#include  <sstream>
#include   <vector>
#include   <time.h>

#include "globals.h"
#include "split_string.h"
#include "precheck_opt.h"

bool ShoreMap_idFilter(int argc, char* argv[])
{
    if(argc < 4)
    {
        printf("\nThis function filters indels from two samples (file format SHOREmap)\n");
        printf("Usage: SHOREmap idFilter [options]                                           \n");
        printf("#Mandatory:                                                                  \n");
        printf("\t--fg-indel     STRING     file of indels to check                          \n");
        printf("\t--bg-indel     STRING     file of indels for removing/keeping              \n");
        //printf("\t--fg-min-cov   INT        minimum coverage on the indels to keep         \n");
        //printf("\t--bg-min-cov   INT        minimum coverage on the indels to keep         \n");        
        printf("#Optional (TODO):                                                            \n");
        printf("\t--indel-size   INT        Size of indels [0] (if 0, cancel output indels). \n");
        printf("\t--min-AF       DOUBLE     minimum AF to record quality reference [0.2]     \n");
        printf("\t-runid         INT        ID of run [1]                                    \n");
        printf("\nExample1: SHOREmap indelFilter --fg-indel ins.A.txt --bg-indel ins.B.txt -runid 4;\n");
        exit(1);
    }
    
    string fg_indel("");
    string bg_indel("");
    long indel_size = 3;
    double freq_min = 0.2;
    int run_id      = 1;
    
    double startT=clock();
    int ic = 2;
    while (ic < argc) 
    {
        if(!strcmp(argv[ic],"--fg-indel"))
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fg_indel        += argv[ic];
            FILE* fp_fg = fopen((char*)fg_indel.c_str(), "r"); 
            if(fp_fg == NULL)
            {
                printf("fg-indel file \'%s\' does NOT exist. Exited.\n", fg_indel.c_str());
                exit(1);
            }
            printf("\tFg-indel file provided:\t\t%s\n", fg_indel.c_str());
            fclose(fp_fg);
        }
        else if(!strcmp(argv[ic],"--bg-indel"))
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            bg_indel        += argv[ic];
            FILE* fp_bg = fopen((char*)bg_indel.c_str(), "r"); 
            if(fp_bg == NULL)
            {
                printf("bg-indel file \'%s\' does NOT exist. Exited.\n", bg_indel.c_str());
                exit(1);
            }
            printf("\tBg-indel file provided:\t\t%s\n", bg_indel.c_str());
            fclose(fp_bg);
        }
        else if(!strcmp(argv[ic],"--indel-size")) // not implemented in real cases yet
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            indel_size = atol(argv[ic]);
            if(indel_size<0) 
                { printf("ERROR: arg of %s must be >= 0. Exited.\n",
                       argv[ic-1]); exit(1);}
            printf("\tindel_size set as:\t\t%ld (TODO)\n", peak_window_step);
        }
        else if(!strcmp(argv[ic],"--min-AF"))    // not implemented in real cases yet
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            freq_min     = atof(argv[ic]);
            if(freq_min<0.0 || freq_min>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            printf("\tfreq_min set as:\t\t\t%.2f. (TODO)\n", reg_freq_min);
        }
        else if(!strcmp(argv[ic],"-runid"))
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            runid = atoi(argv[ic]);
            printf("\trunid set as:\t\t\t\t%ld.\n", runid);
        }
        else
        { 
            printf("\tWarning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        ic ++;
    }   
    
    //reset fg_indel and bg_indel
    string fg_indel_exlpath("");
    string bg_indel_exlpath("");
    vector<string> indelfileinfo = split_string(fg_indel, '/');
    fg_indel_exlpath = indelfileinfo[indelfileinfo.size()-1]; 
    indelfileinfo.clear();
    indelfileinfo = split_string(bg_indel, '/');
    bg_indel_exlpath = indelfileinfo[indelfileinfo.size()-1]; 
    
    // output files: specific indels v.s. common indels
    std::stringstream ss;
    ss.str("");
    ss << runid << "_" << "specific_fg_" << fg_indel_exlpath.c_str() << '\0' << endl;
    fstream fpfgsp;
    fpfgsp.open ((ss.str()).c_str(), ios::out);
    if(!fpfgsp.is_open())
    {
        printf("Cannot open file to write specific indels (in ShoreMap_idFilter.cpp). ");
        printf("Exited. \n");
        exit(1);
    }
    ss.str("");
    ss << runid << "_" << "specific_bg_" << bg_indel_exlpath.c_str() << '\0' << endl;
    fstream fpbgsp;
    fpbgsp.open ((ss.str()).c_str(), ios::out);
    if(!fpbgsp.is_open())
    {
        printf("Cannot open file to write specific indels (in ShoreMap_idFilter.cpp). ");
        printf("Exited. \n");
        exit(1);
    }
    ss.str("");
    ss << runid << "_" << "common_fg_" << fg_indel_exlpath.c_str() << '\0' << endl;
    fstream fpfgcm;
    fpfgcm.open ((ss.str()).c_str(), ios::out);
    if(!fpfgcm.is_open())
    {
        printf("Cannot open file to write common indels (in ShoreMap_idFilter.cpp). ");
        printf("Exited. \n");
        exit(1);
    }
    ss.str("");
    ss << runid << "_" << "common_bg_" << bg_indel_exlpath.c_str() << '\0' << endl;
    fstream fpbgcm;
    fpbgcm.open ((ss.str()).c_str(), ios::out);
    if(!fpbgcm.is_open())
    {
        printf("Cannot open file to write common indels (in ShoreMap_idFilter.cpp). ");
        printf("Exited. \n");
        exit(1);
    }
    // files of indels to read
    std::ifstream fpfg (fg_indel.c_str());
    if(!fpfg.is_open())
    {
        printf("Fg-indel file \'%s\' does NOT exist. Exited.\n", fg_indel.c_str());
        exit(1);
    }
    std::ifstream fpbg (bg_indel.c_str());
    if(!fpbg.is_open())
    {
        printf("Bg-indel file \'%s\' does NOT exist. Exited.\n", bg_indel.c_str());
        exit(1);
    }
    
    map<string, string> bg_indels;
    
    cout << "Start parsing indels from bg file... " << endl;
    while(fpbg.good())
    {
        string line("");
        getline(fpbg, line);
        if(line.size()==0 ||line.find("#")!=std::string::npos) continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        // 0:proj_1	1:chr1	2:28851	3:28853	4:1	5:A	6:UNKWN	7:14	8:0.5185	9:1.0000
        string key = lineinfo[1] + "#" + lineinfo[2] + "#" + lineinfo[5];
        bg_indels.insert(std::pair<string, string>(key, line));
    }
    cout << "Done. Number of indels from bg file " << bg_indels.size() << endl;
    
    cout << "Start parsing indels from fg file... " << endl;
    int ncommon  = 0;
    int nspecfic = 0;
    while(fpfg.good())
    {
        string line("");
        getline(fpfg, line);
        if(line.size()==0 ||line.find("#")!=std::string::npos) continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        string key = lineinfo[1] + "#" + lineinfo[2] + "#" + lineinfo[5];
        map<string, string>::iterator indel_itr;
        indel_itr  = bg_indels.find(key);
        if(indel_itr != bg_indels.end()) // common: but they can have different coverage, AF info etc.
        {
            fpfgcm << line << endl;
            fpbgcm << (*indel_itr).second << endl;
            (*indel_itr).second = "c";   // clear common ones; remainings are specific to bg
            ncommon ++;
        }
        else
        {
            fpfgsp << line << endl;
            nspecfic ++;
        }
    }
    cout << "Indels specific to fg, common between fg and bg have been recorded in " << endl
         << "\t" << runid << "_specific_fg_" << fg_indel_exlpath << "(" << nspecfic << "), " << endl
         << "\t" << runid << "_common_fg_"   << fg_indel_exlpath << "(" << ncommon  << "), " << endl
         << "\t" << runid << "_common_bg_"   << bg_indel_exlpath << "(" << ncommon  << "), " << endl;
    // record specific indels to bg
    int nspecificbg = 0;
    map<string, string>::iterator indel_itr;
    map<string, string>::iterator indel_itr_end;
    indel_itr     = bg_indels.begin();
    indel_itr_end = bg_indels.end();
    while(indel_itr != indel_itr_end)
    {
        if((*indel_itr).second != "c") 
        {
            fpbgsp << (*indel_itr).second << endl;
            nspecificbg ++;
        }
        indel_itr ++;
    }
    cout << "\t" << runid  << "_specific_bg_"  << bg_indel_exlpath       << "(" << nspecificbg << "), "
         << endl           << endl             << "respectively."        << endl;
    
    fpfgsp.close();
    fpbgsp.close();
    fpfgcm.close();
    fpbgcm.close();
    fpfg.close();
    fpbg.close();
    
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nidFilter function successfully finished on %s\n", asctime(tinfo));
    return true;
}
