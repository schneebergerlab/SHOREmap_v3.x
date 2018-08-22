/* this function creates markers for F2 by using SNPs of parent(s). Those parental SNPs, which 
indicate marker-positions, are checked by using background ref-base information. Then all the SNPs
are combined with SNPs of F2 to obtain markers. Assuming F2 and parent a have the same phenotype
while parent b does not have the feature:

   basic scheme: 1.check SNP of parent a by using ref-base information of parent b;
                 2.check SNP of parent b by using ref-base information of parent a (for this parent,
                       a SNP can also be chekced by using ref-base information of F2???)
                 3.
   TODO: how about directly obtaining allele frequency of marker-positions with consensus 
         information of F2?                               
*/
#include  <stddef.h>
#include  <stdlib.h>
#include   <stdio.h>
#include    <math.h>
#include    <string>
#include  <string.h>
#include   <ctype.h>
#include    <time.h>
#include       <map>
#include  <iostream>

// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>

#include "globals.h"

#include "is_number.h"
#include "precheck_opt.h"
#include "print_error_exit.h"
#include "read_chromosomes.h"
//#include "find_marker_pos_2parents_v3.h" // 2013-09-17
#include "find_marker_pos_2parents_v3.h"   // 2014-03-28

int  cmd_create(int argc, char* argv[]);

bool ShoreMap_create2(int argc, char* argv[])
{
    if (argc < 10) 
    {
        std::string create_usage = "";
        printf("Usage: SHOREmap create [Options]\n\n");
        printf("#Mandatory:\n");
        printf("--chrsizes              STRING   File of name and size of chromosomes;   [NULL]\n");
        printf("--folder                STRING   Output folder (with path);              [NULL]\n");
        printf("--marker                STRING   SNP-file of targeted mutant (F2);       [NULL]\n");
        printf("--marker-pa             STRING   SNP-file of bg (parent) a;              [NULL]\n");
        printf("                                 with the same phenotype as the mutant;        \n");
        printf("--marker-pb             STRING   SNP-file of bg (parent) b;              [NULL]\n");
        printf("                                 without the same phenotype as the mutant;     \n");
        printf("--bg-ref-base-pa        STRING   Base-call file of bg a;                 [NULL]\n");
        printf("                                 with   the same phenotype as the mutant.      \n");        
        printf("--bg-ref-base-pb        STRING   Base-call file of bg b;                 [NULL]\n");
        printf("                                 without the same phenotype as the mutant.   \n\n");
        printf("#Optional:\n");     
        /*quality, coverage, allele frequency of fg: targeted SNPs                                */
        //'N' and '-' information are not used yet!!! If SNPs of parent x are to be filtered, then
        //then the consensus information of px-SNPs must be provided to do this filtering.
        //printf("--fg-N-cov            INT      Maximum coverage    of foreground 'N'-call;\n"); 
        //printf("                              (default: INF).\n");	
        //printf("--fg-INDEL-cov        INT      Maximum coverage    of foreground indel-call;\n"); 
        //printf("                               (default: INF).\n");
        printf("#thresholds on background (parent a/b) SNPs for creating markers in the pool:\n");
        printf("--pmarker-min-cov       INT      Minimum coverage  of bg SNPs               [0]\n");
        printf("--pmarker-max-cov       INT      Maximum coverage  of bg SNPs             [INF]\n");
        printf("--pmarker-min-freq      DOUBLE   Minimum frequency of bg SNPs             [0.0]\n");
        printf("--pmarker-score         INT      Minimum quality   of base-call of bg SNPs  [0]\n");
        printf("#constraints on background (parent b/a) ref-base-call to filter (a/b) SNPs:\n");
        printf("--bg-ref-cov            INT      Minimum coverage  of bg ref base call;    [10]\n");
        printf("--bg-ref-cov-max        INT      Maximum coverage  of bg ref base call;   [500]\n");
        printf("--bg-ref-freq           DOUBLE   Minimum frequency of bg ref base call;   [0.8]\n");
        printf("--bg-ref-score          INT      Minimum quality   of bg ref base call;    [25]\n");
        printf("--keep-common           Keep common snps of parental lines in markers list[off]\n");
        printf("-verbose                         Be talkative.                            [off]\n");
        printf("*bg-background.                                                                \n");        
        exit(1);
    }
    
    /* set defaults */
    /*  defaults */
    fchrsizes           =  "";   // CMD["--chrsizes"]
    out_folder          =  "";   // CMD["--folder"]
    fmarker             =  "";   // CMD["--marker"]
    fmarker_pa          =  "";   // CMD["--marker-pa"]
    fmarker_pb          =  "";   // CMD["--marker-pb"]
    bg_ref_base_file    =  "";
    bg_ref_base_file_pb =  "";
    pmarker_score       =   0;   // CMD["--pmarker-score"]
    pmarker_min_cov     =   0;   // CMD["--pmarker-min-cov"]
    pmarker_max_cov     = INF;   // CMD["--pmarmer-max-cov"]
    pmarker_min_freq    = 0.0;   // CMD["--pmarker-min-freq"]
    bg_ref_cov          =  10;
    bg_ref_cov_max      = 500;
    bg_ref_freq         = 0.8;
    bg_ref_score        =  25;
    fg_N_cov            = INF;   // <=> marker_N_cov     - not used;// INF means no filtering
    fg_INDEL_cov        = INF;   // <=> marker_INDEL_cov - not used;// INF means no filtering
    filter_min_coverage =   0;   // will be set by pmarker_min_cov  // 0   means no filtering
    filter_max_coverage = INF;   // will be set by pmarker_max_cov  // INF means no filtering
    reg_freq_min        =   0;   // will be set by pmarker_min_freq // 0   means no filtering
    marker_score        =   0;   // will be set by pmarker_score    // 0   means no filtering  
    keep_common         = false;
    verbose             =   0;   // CMD["-verbose"]
    /* read command options */
        
    cmd_create(argc, argv);
    /* the following parameters are used in check_ref_err.cpp - to be consistent with backcross   */
    filter_min_coverage = pmarker_min_cov;
    filter_max_coverage = pmarker_max_cov;
    reg_freq_min        = pmarker_min_freq;
    marker_score        = pmarker_score;
    std::string  F2marker_ParentAsumBminusAB;
    
    /* find_marker_pos_2parents_v2 only uses resequencing data of parental lines - 2013-09-17 10:01*/
    /* find_marker_pos_2parents_v3 uses resequencing data of parental lines&F2   - 2014-03-28 19:38*/
    find_marker_pos_2parents_v3( (char*) fmarker_pa.c_str(), 
                                 (char*) bg_ref_base_file_pb.c_str(), 
                                 (char*) fmarker_pb.c_str(),
                                 (char*) bg_ref_base_file.c_str(),
                                 (char*) fmarker.c_str(),
                                 bg_ref_cov,
                                 bg_ref_cov_max,
                                 bg_ref_freq,             /* minimum freq */
                                 fg_N_cov,
                                 fg_INDEL_cov, 
                                 bg_ref_score,
                                 &F2marker_ParentAsumBminusAB);
    cout << "Created markers for analysing the pool has been recorded: \t";
    cout << F2marker_ParentAsumBminusAB << endl;
}

int  cmd_create(int argc, char* argv[])
{
    /* 
       Read practical values from cmd options; option-values are checked and updated
       in map<string, string> CMD with format: CMD["option"] = "arg" (if exist).
    */
    int ic = 2;
    while (ic < argc)
    {
        if(!strcmp(argv[ic],"-verbose"))      // option variable: verbose
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            verbose         = 1;
            printf("Be talkative during process. \n");
        }
        ic ++;
    }   
    ic = 2;                                   // option ic=0: SHOREmap; ic=1: extract
    while (ic < argc) 
    {
        if(!strcmp(argv[ic],"--chrsizes"))    // option 1 to variable: fchrsizes
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fchrsizes        += argv[ic];     // set value for option
            FILE* fp_chrsizes = fopen((char*)fchrsizes.c_str(), "r");
            if(!fp_chrsizes)  { print_error_exit(argv[ic-1], false); }
            fclose(fp_chrsizes);              // check done.
            if (verbose) printf("Chrs sizes read from file:\t\t%s\n", fchrsizes.c_str());
            if(!read_chromosomes((char*)fchrsizes.c_str()))  // read contents
            {
                printf("ERROR: invalid content in file %s.\n", fchrsizes.c_str());
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--folder")) // option 2 to variable: out_folder
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            out_folder         += argv[ic];
            if(out_folder[out_folder.length()-1] != '/') out_folder += "/";
            DIR* dir_out_folder = opendir((char*)out_folder.c_str());
            if(dir_out_folder  == NULL)
            {
                if(!mkdir((char*)out_folder.c_str(), S_IRWXU|S_IRWXG|S_IRWXO))
                {
                    /* if !mkdir() is TRUE: a new directory has to be created.
                    if out_folder is "outfolder/" instead of "/your/output/path/outfolder", 
                    an "outfolder/" will be created under the current working directory */
                    if (verbose) printf("Folder created:\t\t\t\t\t%s.\n", out_folder.c_str());
                }
                else 
                {
                    printf("ERROR: cannot create output path. Exited.\n");
                    exit(1);
                }
            }
            closedir(dir_out_folder);
        }
        else if(!strcmp(argv[ic],"--marker"))  // option 3 to variable: fmarker
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker        += argv[ic];
            FILE* fp_marker = fopen((char*)fmarker.c_str(), "r"); 
            if(fp_marker == NULL)
            {
                printf("marker file \'%s\' does NOT exist. Exited.\n", fmarker.c_str());
                exit(1);
            }
            if (verbose) printf("File of markers provided:\t\t%s\n", fmarker.c_str());
           // read_marker();
            fclose(fp_marker);
        }
        else if(!strcmp(argv[ic],"--marker-pa"))  // option 4.1 to variable: fmarker_pa
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker_pa        += argv[ic];
            FILE* fp_marker_pa = fopen((char*)fmarker_pa.c_str(), "r"); 
            if(fp_marker_pa == NULL)
            {
                printf("marker file of parent a \'%s\' does NOT exist. Exited.\n", fmarker_pa.c_str());
                exit(1);
            }
            if (verbose) printf("File of markers provided as parent a:\t%s\n", fmarker_pa.c_str());
            fclose(fp_marker_pa);
        }
        else if(!strcmp(argv[ic],"--marker-pb"))  // option 4.2 to variable: fmarker_pb
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker_pb        += argv[ic];
            FILE* fp_marker_pb = fopen((char*)fmarker_pb.c_str(), "r"); 
            if(fp_marker_pb == NULL)
            {
                printf("marker file of parent b \'%s\' does NOT exist. Exited.\n", fmarker_pb.c_str());
                exit(1);
            }
            if (verbose) printf("File of markers provided as parent b:\t%s\n", fmarker_pb.c_str());
            fclose(fp_marker_pb);
        }
        else if(!strcmp(argv[ic],"--bg-ref-base-pa"))  // option 4.3 to variable: bg_ref_base_file
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            bg_ref_base_file = argv[ic];
            FILE* fp_bgrbf   = fopen((char*)bg_ref_base_file.c_str(), "r"); 
            if(fp_bgrbf     == NULL)
            {
                printf("consensus file \'%s\' does NOT exist. Exited.\n", bg_ref_base_file.c_str());
                exit(1);
            }
            fclose(fp_bgrbf);
        }
        else if(!strcmp(argv[ic],"--bg-ref-base-pb"))  // option 4.4 to variable: bg_ref_base_file_pb
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            bg_ref_base_file_pb = argv[ic];
            FILE* fp_bgrbf   = fopen((char*)bg_ref_base_file_pb.c_str(), "r"); 
            if(fp_bgrbf     == NULL)
            {
                printf("consensus file \'%s\' does NOT exist. Exited.\n", bg_ref_base_file_pb.c_str());
                exit(1);
            }
            fclose(fp_bgrbf);
        }
        else if(!strcmp(argv[ic],"--bg-ref-cov"))     // option 4.5 to variable: bg_ref_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov set as:\t\t\t%ld\n", bg_ref_cov);
        }
        else if(!strcmp(argv[ic],"--bg-ref-cov-max")) // option 4.6 to variable: bg_ref_cov_max
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov_max   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov_max set as:\t\t\t%ld\n", bg_ref_cov_max);
        }
        else if(!strcmp(argv[ic],"--bg-ref-freq"))   // option 4.7 to variable: bg_ref_freq
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_freq     = atof(argv[ic]);
            if(bg_ref_freq<0.0 || bg_ref_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("bg_ref_freq:\t\t\t\t%.2f.\n", bg_ref_freq);
        }
        else if(!strcmp(argv[ic],"--bg-ref-score"))  // option 4.8 to variable: bg_ref_score
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_score           = atol(argv[ic]);
            if(bg_ref_score<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("bg_ref_score:\t\t\t\t%ld\n", bg_ref_score);
        }
        else if(!strcmp(argv[ic],"--pmarker-score"))  // option 5.1 to variable: pmarker_score
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_score           = atol(argv[ic]);
            if(pmarker_score < 1)
                { printf("ERROR: arg of %s must be >= 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("pmarker_score:\t\t\t\t%ld\n", pmarker_score);
        }
        else if(!strcmp(argv[ic],"--pmarker-min-cov"))// option 5.2 to variable: pmarker_min_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_min_cov   = atol(argv[ic]);
            if (verbose) printf("pmarker_min_cov set as:\t\t\t%ld\n", pmarker_min_cov);
        }
        else if(!strcmp(argv[ic],"--pmarker-max-cov"))// option 5.3 to variable: pmarker_max_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_max_cov   = atol(argv[ic]);
            if (verbose) printf("pmarker_max_cov set as:\t\t\t%ld\n", pmarker_max_cov);
        }
        else if(!strcmp(argv[ic],"--pmarker-min-freq")) // option 5.4 to variable: pmarker_min_freq
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_min_freq     = atof(argv[ic]);
            if(pmarker_min_freq<0.0 || pmarker_min_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            if (verbose) printf("pmarker_min_freq set as:\t\t%.2f.\n", pmarker_min_freq);
        }
        else if(!strcmp(argv[ic],"--keep-common"))     // option 5.5 to variable: keep_common
        {
            keep_common = true;
            if (verbose) printf("keep_common set as:\t\ttrue.\n");
        }
        else // other options not recognized
        { 
             if(!strcmp(argv[ic],"-verbose"));
             else printf("Warning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        /* check mandatory options */
        if(ic == argc-1)
        {
            if(CMD.find("--chrsizes")==CMD.end())
            { 
                printf("ERROR: chr sizes file           is NOT provided. Exited.\n"); exit(1);
            }
            if (CMD.find("--folder")==CMD.end()) 
            {
                printf("ERROR: output folder            is NOT provided. Exited.\n"); exit(1);
            }
            if (CMD.find("--marker-pa")==CMD.end() && CMD.find("--marker-pb")==CMD.end())
            {
                printf("ERROR: --marker-pa and/or --marker-pb must be provided. Exited. \n"); 
                exit(1);
            } 
            if(CMD.find("--marker-pa")!=CMD.end() && CMD.find("--bg-ref-base-pb")==CMD.end())  
            {
                printf("ERROR: --marker-pa is provided, requiring --bg-ref-base-pb. Exited.\n");
                exit(1);                
            }   
            if(CMD.find("--marker-pb")!=CMD.end() && CMD.find("--bg-ref-base-pa")==CMD.end())  
            {
                printf("ERROR: --marker-pb is provided, requiring --bg-ref-base-pa. Exited. \n");
                exit(1);                
            } 
            if(CMD.find("--marker") == CMD.end())
            {
                printf("WARNING: --marker missed! Clustering will not function properly!\n");
            }
        }
        ic ++;
    }
    return 1;
}
