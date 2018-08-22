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
#include "read_chromosomes.h"
#include "read_marker.h"
#include "read_allele_count2.h"
#include "precheck_opt.h"
#include "print_error_exit.h"
#include "read_extract_cons.h"
#include "read_extract_bg_qref.h"


int  cmd_extractC(int argc, char* argv[]);

int ShoreMap_extract(int argc, char* argv[])
{
     /* print out usage of outcross */
    if (argc < 4) {
        printf("\nUsage: ./SHOREmap extract [Options].\n\n");
        printf("Mandatory:\n");
        printf("--chrsizes              STRING   File of chromosome names and sizes.\n");
        printf("--folder                STRING   Output folder.\n");
        printf("--marker                STRING   File of candidate markers\n");
        printf("--consen                STRING   Consensus file.\n");
        printf("--extract-bg-ref                 on : from bg-quality-reference file.\n");
        printf("                                 off: from fg-consensus-summary file(default).\n");
        //printf("Optionals (// TODO - 2013-03-23 15:06):\n");
        //printf("row/line and column interested:\n"); 
        //printf("--row-first             INT      start of line of consensus.\n");
        //printf("                                 (default: 1).\n");
        //printf("--row-last              INT      end of line of consensus.\n");
        //printf("                                 (default: 0: END_OF_FILE).\n");
        //printf("--col-first             INT      start of column of consensus.\n");
        //printf("                                 (default: 1).\n");
        //printf("--col-last              INT      end of column of consensus.\n");
        //printf("                                 (default: 8).\n");
        printf("-verbose                         Be talkative.\n\n");
        exit(1);
    }
    
    /* initial default */
    fchrsizes  = "";      //CMD["--chrsizes"]    = "";
    out_folder = "";      //CMD["--folder"]      = "";
    fmarker    = "";      //CMD["--marker"]      = "";
    fconsensus = "";      //CMD["--consen"]      = "";
    extract_bg_ref = 0;
    row_first  = 0;       // 0: parse from the first line of a file
    row_last   = 0;       // 0: parse until the last line of a file
    col_first  = 1;       // 1.chr       2.position  3.cons_base 4.coverage  (if fg-consen)
    col_last   = 8;       // 5.A_support 6.C_support 7.G_support 8.T_support (if fg-consen)
    marker_score = 0;
    
    /* read parameter settings from cmd line */
    cmd_extractC(argc, argv);

    /* read markers */
    unsigned long num_m_given = 0;
    unsigned long num_m_inuse = 0;
    if(!read_marker((char*)fmarker.c_str(), &num_m_given, &num_m_inuse))
    {
       printf("ERROR: no marker info recorded. Exited. \n"); 
       printf("Hint 1: check if chromosome ids of all files are consistent.\n");
       printf("Hint 2: check if files of candidate markers and consensus calls are consistent.\n");
       exit(1);
    }
    /* read and extract fg-consensus/bg-ref-base info */
    bool read_flag = false; // if successfully read at least 1 marker, it is reset as true
    if(extract_bg_ref == 1)
    {
        // e.g., from quality_reference.txt of bg-mutant
        read_flag = read_extract_bg_qref((char*)fconsensus.c_str());
    }
    else
    {
        // e.g., from consensus_summary.txt of fg-mutant
        read_flag = read_extract_cons((char*)fconsensus.c_str());
    }
    /* write cmd log */
    std::string fETlog("");
    fETlog = out_folder +  "Extracting.log";
    FILE* fplog = fopen(fETlog.c_str(), "a+");
    if(fplog)
    {
        fprintf(fplog, "Inputs provided:\n\n");
        map<string, string>::iterator cmditr = CMD.begin();
        while(cmditr != CMD.end())
        {
            if((*cmditr).first.length()<8)
            {
                fprintf(fplog, "%s \t\t%s\n", (*cmditr).first.c_str(), (*cmditr).second.c_str());
            }
            else
            {
                fprintf(fplog, "%s \t%s\n", (*cmditr).first.c_str(), (*cmditr).second.c_str());
            }
            cmditr++;
        }
    }
    else
    {
        printf("Cannot open file to write log. Exited (in init_backcross).\n"); exit(1);
    }
    if(read_flag)
    {
        fprintf(fplog, "\n\nOutputs achieved:\t");
        if(extract_bg_ref == 0)
        fprintf(fplog, "%sextracted_consensus_%ld.txt\n", out_folder.c_str(), row_first);
        else
        fprintf(fplog, "%sextracted_quality_ref_base_%ld.txt\0", out_folder.c_str(), row_first);
    }
    else
    {
        fprintf(fplog, "\n\nWarning: no consensus info recorded for given markers. \n");
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    fprintf(fplog, "\nExtract function successfully finished on %s\n", asctime(tinfo));
    fclose(fplog);
    
    return 0;
}

int  cmd_extractC(int argc, char* argv[])
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
            fclose(fp_marker);
        }
        else if(!strcmp(argv[ic],"--consen"))  // option 4 to variable: fconsensus
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fconsensus        += argv[ic];
            FILE* fp_consensus = fopen((char*)fconsensus.c_str(), "r"); 
            if(fp_consensus   == NULL)
            {
                printf("consensus file \'%s\' does NOT exist. Exited.\n", fconsensus.c_str());
                exit(1);
            }
            fclose(fp_consensus);
        }
        else if(!strcmp(argv[ic],"--row-first"))  // option 5 to variable: row_first
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            row_first   = atol(argv[ic]);
            if (verbose) printf("row_first set as:\t\t%ld\n", row_first);
        }
        else if(!strcmp(argv[ic],"--row-last"))  // option 6 to variable: row_last
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            row_last   = atol(argv[ic]);
            if (verbose) printf("row_last set as:\t\t%ld\n", row_last);
        }
        else if(!strcmp(argv[ic],"--col-first")) // option 7 to variable: col_first
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            col_first   = atol(argv[ic]);
            if (verbose) printf("col_first set as:\t\t%ld\n", col_first);
        }
        else if(!strcmp(argv[ic],"--col-last"))  // option 8 to variable: col_last
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            col_last   = atol(argv[ic]);
            if (verbose) printf("col_last set as:\t\t%ld\n", col_last);
        }
        else if(!strcmp(argv[ic],"--extract-bg-ref"))  // option 9 to variable: extract_bg_ref
        {
            extract_bg_ref = 1;                        
            if (verbose) printf("extract_bg_ref set as:\t\t%ld\n", extract_bg_ref);
        }
        else // other options not recognized
        { 
             if(!strcmp(argv[ic],"-verbose"));
             else printf("Warning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        /* check mandatory options */
        if(ic == argc-1)
        {
            if(CMD.find("--chrsizes")==CMD.end()){ 
                printf("ERROR: chr sizes file           is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--folder")==CMD.end()) {
                printf("ERROR: output folder            is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--marker")==CMD.end()) {
                printf("ERROR: marker file              is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--consen")==CMD.end()) {
                printf("ERROR: fg-consensus/bg-ref file is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--extract-bg-ref") == CMD.end()){
                printf("Extract from a consensus file.\n");                                   }
        }
        ic ++;
    }
    return 1;
}
