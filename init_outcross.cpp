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
#include "read_referror.h"
#include "read_marker.h"
#include "read_allele_count.h"
#include "read_allele_count2.h"
#include "write_log.h"
#include "write_zoom_region.h"
#include "precheck_opt.h"
#include "print_error_exit.h"

//declare functions
int  cmd_outcross(int argc, char* argv[]);

/* this function initialize the outcross */
int init_outcross(int argc, char* argv[])
{
    /* print out usage of outcross */
    if (argc < 10) {
        printf("Usage: SHOREmap outcross [Options].\n\n");
        printf("Mandatory:\n");
        printf("--chrsizes              STRING   Tabed file with chromosome name and size.\n");
        printf("--folder                STRING   Output folder.\n");
        printf("--marker                STRING   Marker file.\n");
        printf("--consen                STRING   Consensus file.\n\n");
        printf("Optionals:\n");
        printf("Confidence interval:.\n");//
        printf("-conf-int                        Switch on confidence interval calculation.\n");
        printf("                                 (default: 1.0).\n");
        printf("--conf                  DOUBLE   Confidence level.\n");
        printf("                                 (default: 0.99).\n");
        printf("--peak-window-size      INT      Size of initial window.\n");
        printf("                                 (default: 50000).\n");
        printf("--peak-window-step      INT      Distance of initial windows.\n");
        printf("                                 (default: 10000).\n");
        printf("--mis-phenotyped        DOUBLE   Degree of putatively mis-scored plants.\n");
        printf("                                 (default: 0.00).\n\n");
        printf("Visulization:\n");     //
        printf("--window-size           INT      (default: 50000).\n");
        printf("--window-step           INT      Used for smoothed visulization.\n");
        printf("                                 and \"boost\"-value calculation.\n");
        printf("                                 (default: 10000).\n");
        printf("-plot-marker                     Plot single markers.\n");
        printf("-plot-r                          Plot frequency calculation (\"r\");).\n");
        printf("                                 instead of \"boost\".\n\n");
        printf("Filter:\n");           //
        printf("--min-marker            INT      Filter windows with low numbers.\n");
        printf("                                 of markers (default: 10).\n");
        printf("--min-coverage          INT      Filter single marker with low.\n");
        printf("                                 average coverage (default: 0).\n");
        printf("--max-coverage          INT      Filter single marker with high.\n");
        printf("                                 average coverage (default: Inf).\n");
        printf("-outlier                         Switch on outlier filtering.\n");
        printf("--outlier-window-size   INT      Window size to assess local.\n");
        printf("                                 allele frequency used for.\n");
        printf("                                 outlier removal.\n");
        printf("                                 (default: 200000).\n\n");
        printf("--outlier-pvalue        DOUBLE   p-value used for outlier.\n");
        printf("                                 removal (default: 0.05).\n\n");
        printf("Zooming:\n");          //
        printf("--chromosome            INT      Zoom to chromosome ...\n");
        printf("--begin                 INT      .. from here ...\n");
        printf("--end                   INT      .. to here with a ...\n");
        printf("--minfreq               INT      .. minimal to ...\n");
        printf("--maxfreq               INT      .. maximal frequency..\n\n");
        printf("--referrors             STRING   Reference errors file.\n");
        printf("-background2                     Mutation is in second parent.\n\n");
        printf("-verbose                         Be talkative.\n");
        printf("See documentation for file formats..\n");
        exit(1);
    }
    
    /* initial default */
    expect              = 1.0;     //CMD["-target"]       = "";
    confidence          = 2;       //CMD["-conf-int"]     = ""; 
                                   //CMD["--conf"]        = "";
    fchrsizes           = "";      //CMD["--chrsizes"]    = "";
    out_folder          = "";      //CMD["--folder"]      = "";
    fmarker             = "";      //CMD["--marker"]      = "";
    fconsensus          = "";      //CMD["--consen"]      = "";
    window_size         = 50000;   //CMD["--window-size"] = "";
    window_step         = 10000;   //CMD["--window-step"] = "";

    peak_window_size    = 50000;   //CMD["--peak-window-size"] = "";
    peak_window_step    = 10000;   //CMD["--peak-window-step"] = "";

    filter_min_marker   = 10;      //CMD["--min-marker"]  = "";
    filter_min_coverage = 10;      //CMD["--min-coverage"]= ""; // caution: case-sensitive
    filter_max_coverage = INF;     //CMD["--max-coverage"]= "";

    outlier_window_size = 200000;  //CMD["--outlier-window-size"] = "";
    outlier_pvalue      = 0.05;    //CMD["--outlier-pvalue"]      = "";
    misphenotyped       = 0.00;    //CMD["--mis-phenotyped"]      = "";

    reg_chromosome      = "";      //CMD["--chromosome"] = "";
    reg_begin           = 0;       //CMD["--begin"]      = "";
    reg_end             = 0;       //CMD["--end"]        = "";
    reg_freq_min        = 0;       //CMD["--minfreq"]    = "";
    reg_freq_max        = 0;       //CMD["--maxfreq"]    = "";

    freferror           = "";      //CMD["--referrors"]  = "";
    background2         = 0;       //CMD["-background2"] = "";
    verbose             = 0;       //CMD["-verbose"]     = "";
    
    boost_max           = 10000;   //CMD["-boost-max"]   = "";
    r_max               = 10000;   //CMD["-r-max"]       = "";
    plot_boost          = 1;       // related to "-plot-r"
    plot_r              = 0;       //CMD["-plot-r"]      = "";
    plot_marker         = 0;       //CMD["-plot-marker-off"] = "";
                             
    runid               = 1;       //CMD["-runid"]       = "";

    /* 0.read practical values from cmd line */
    cmd_outcross(argc, argv);

    /* 1.read reference error positions (optional) */
    if(CMD.find("--referrors") != CMD.end())
    {
        read_referror((char*)freferror.c_str());
    }
    /* 2.read markers */
    unsigned long num_m_given = 0;
    unsigned long num_m_inuse = 0;
    if(!read_marker((char*)fmarker.c_str(), &num_m_given, &num_m_inuse))
    {
       printf("ERROR: no marker info recorded. Exited. \n"); exit(1);
    }
    
    /* test another reading */
    //if (!read_allele_counts2((char*)fconsensus.c_str()))
    //{
    //    printf("ERROR: no consensus info recorded. Exited. \n"); exit(1);
    //    exit(1);
    //}
    
    /* end of test */
    
    /* 3. read consensus base and allele1 allele2 error counts */
    if (!read_allele_counts2((char*)fconsensus.c_str()))
    {
        printf("ERROR: no consensus info recorded. Exited. \n"); exit(1);
    }

    return 1;
}

////////////////////////////////////////
int cmd_outcross(int argc, char* argv[])
{
    /* 
       Read practical values from cmd options; option-values are checked and updated
       in map<sstring, string> CMD with format: CMD["option"] = "arg" (if exist).
    */
    int ic = 2; // check verbose first.
    while (ic < argc)
    {
        if(!strcmp(argv[ic],"-verbose"))      // option 27 to variable: verbose
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            verbose         = 1;
            printf("Be talkative during process. \n");
        }
        ic ++;
    }
    ic = 2;    // option ic=0: SHOREmap; ic=1: outcross
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
            DIR* dir_out_folder = opendir((char*)out_folder.c_str());
            if(dir_out_folder  == NULL)
            {
                if(!mkdir((char*)out_folder.c_str(), S_IRWXU|S_IRWXG|S_IRWXO))
                {
                 /* if !mkdir() is TRUE: a new directory has to be created.
                    if out_folder is "outfolder/" instead of "/your/output/path/outfolder", 
                    an "outfolder/" will be created under the current working directory */
                    if (verbose) printf("Folder\t\t\t\t\t%s created.\n", out_folder.c_str());
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
            // read_allele_counts();
            fclose(fp_consensus);
        }
        else if(!strcmp(argv[ic],"-conf-int"))          // option 5 to variable: confidence
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            confidence       = atof(argv[ic]);
            if (verbose) printf("Confidence set as:\t\t\t%.2f by option %s.\n", confidence, argv[ic-1]);
        }
        else if(!strcmp(argv[ic],"--conf"))             // option 6 to variable: confidence
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            confidence    = atof(argv[ic]);
            if(confidence<0.0 || confidence>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            if (verbose) printf("confidence set as:\t\t\t%.2f by option %s (reset: %d).\n", 
            confidence, argv[ic-1], CMD.find("-conf-int")!=CMD.end());
        }
        else if(!strcmp(argv[ic],"--peak-window-size")) // option 7 to variable: peak_window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            peak_window_size          = atol(argv[ic]);
            if(peak_window_size<1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("peak_window_size set as:\t\t%ld\n", peak_window_size);
        }
        else if(!strcmp(argv[ic],"--peak-window-step")) // option 8 to variable: peak_window_step
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            peak_window_step          = atol(argv[ic]);
            if(peak_window_step<1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("peak_window_step set as:\t\t%ld\n", peak_window_step);
        }
        else if(!strcmp(argv[ic],"--mis-phenotyped"))  // option 9 to variable: misphenotyped
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            misphenotyped           = atof(argv[ic]);
            if(misphenotyped<0.0 || misphenotyped>1.0) 
                { printf("ERROR: arg of %s must be [0.0, 1.0].\n", argv[ic-1]); exit(1);}
            if (verbose) printf("misphenotyped set as:\t\t\t%.2f\n", misphenotyped);
        }
        else if(!strcmp(argv[ic],"--window-size"))    // option 10 to variable: window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_size           = atol(argv[ic]);
            if(window_size<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("window_size set as:\t\t\t%ld\n", window_size);
        }
        else if(!strcmp(argv[ic],"--window-step"))   // option 11 to variable: window_step
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_step           = atol(argv[ic]);
            if(window_step<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("window_step set as:\t\t\t%ld\n", window_step);
        }
        else if(!strcmp(argv[ic],"-plot-marker-off"))    // option 12 to variable: plot_marker
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_marker         = 0;
            if (verbose) printf("plot_marker set as:\t\t\t%d\n", plot_marker);
        }
        else if(!strcmp(argv[ic],"-plot-r"))         // option 13 to variable: plot_r
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_r         = 1;
            plot_boost     = 0; // caution
            if (verbose) printf("plot_r, plot_boost set as:\t\t%d, %d\n", plot_r, plot_boost);
        }
        else if(!strcmp(argv[ic],"--min-marker"))    // option 14 to variable: filter_min_marker
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_marker   = atol(argv[ic]);
            if(filter_min_marker<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("filter_min_marker set as:\t\t%ld\n", filter_min_marker);
        }
        else if(!strcmp(argv[ic],"--min-coverage"))  // option 15 to variable: filter_min_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_min_coverage set as:\t\t%ld\n", filter_min_coverage);
        }
        else if(!strcmp(argv[ic],"--max-coverage"))  // option 16 to variable: filter_max_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_max_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_max_coverage set as:\t\t%ld\n", filter_max_coverage);
        }
        else if(!strcmp(argv[ic],"-outlier"))       // option 17 to variable: outlier_window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            outlier_window_size = 0;
            if (verbose) printf("outlier removal is on. \n");
        }
        else if(!strcmp(argv[ic],"--outlier-window-size")) // option 18 to variable: outlier_window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            outlier_window_size          = atol(argv[ic]);
            if (verbose) printf("outlier_window_size set as:\t\t%ld\n", outlier_window_size);
        }
        else if(!strcmp(argv[ic],"--outlier-pvalue"))      // option 19 to variable: outlier_pvalue
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            outlier_pvalue          = atof(argv[ic]);
            if(outlier_pvalue<0.0 || outlier_pvalue>1.0)
                { printf("ERROR: arg of %s must be set in [0.0, 1.0]. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("outlier_pvalue set as:\t\t\t%.2f.\n", outlier_pvalue);
        }
        else if(!strcmp(argv[ic], "--chromosome"))        // option 20 to variable: reg_chromosome
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            reg_chromosome     += argv[ic];
            if (verbose) printf("reg_chromosome id set as:\t\t%s.\n", reg_chromosome.c_str());
            if(CHR2SIZE.find(reg_chromosome)==CHR2SIZE.end()) 
            {
                printf("ERROR: chromosome %s to zoom is not found in the list. Exited.\n", 
                    (char*)reg_chromosome.c_str()); 
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--begin"))             // option 21 to variable: reg_begin
        {
            if(reg_chromosome == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_begin      = atol(argv[ic]);
            string schr(reg_chromosome);
            if(reg_begin<0 || reg_begin>CHR2SIZE[schr])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                     argv[ic-1]); exit(1);}
            if (verbose) printf("reg_begin set as:\t\t\t%ld\n", reg_begin);
        }
        else if(!strcmp(argv[ic],"--end"))               // option 22 to variable: reg_end
        {
            if(reg_chromosome == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_end      = atol(argv[ic]);
            if(reg_end<0 || reg_end>CHR2SIZE[reg_chromosome])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("reg_end set as:\t\t\t\t%ld\n", reg_end);
        }
        else if(!strcmp(argv[ic],"--minfreq"))          // option 23 to variable: reg_freq_min
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_freq_min     = atof(argv[ic]);
            if(reg_freq_min<0.0 || reg_freq_min>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            if (verbose) printf("reg_freq_min set as:\t\t\t%.2f.\n", reg_freq_min);
        }
        else if(!strcmp(argv[ic],"--maxfreq"))         // option 24 to variable: reg_freq_max
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_freq_max     = atof(argv[ic]);
            if(reg_freq_max<0.0 || reg_freq_max>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("reg_freq_max set as:\t\t\t%.2f.\n", reg_freq_max);
        }
        else if(!strcmp(argv[ic],"--referrors"))       // option 25 to variable: freferror
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            freferror         += argv[ic];
            FILE* fp_referror  = fopen((char*)freferror.c_str(), "r"); 
            if(fp_referror == NULL)
                { printf("Ref-error file \'%s\' does NOT exist. Exited.\n", freferror.c_str()); exit(1);}
            fclose(fp_referror);
        }
        else if(!strcmp(argv[ic],"-background2"))      // option 26 to variable: background2
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            background2 = 1;
            if (verbose) printf("mutation in second parent has been told. \n");
        }
        else if(!strcmp(argv[ic],"-runid"))            // option 27 to variable: runid
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            runid = atoi(argv[ic]);
            if (verbose) printf("runid has set as:\t\t\t%ld.\n", runid);
        }
        else if(!strcmp(argv[ic],"-boost-max"))       // option 28 to variable: boost_max
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            boost_max         = atof(argv[ic]);
            if (verbose) printf("boost_max has set as:\t\t\t%.2f.\n", boost_max);
        }
        else if(!strcmp(argv[ic],"-r-max"))           // option 29 to variable: r_max
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            r_max         = atof(argv[ic]);
            if (verbose) printf("r_max has set as:\t\t\t%ld.\n", r_max);
        }
        else if(!strcmp(argv[ic],"-target"))           // option 30 to variable: expect
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            expect         = atof(argv[ic]);
            if (verbose) printf("target/expect has set as:\t\t%.2f.\n", expect);
        }
        else // other options not necessary
        { 
             if(!strcmp(argv[ic],"-verbose"));         // this has been set
             else printf("Warning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        
        if(ic == argc-1) // check mandatory options
        {
            if(CMD.find("--chrsizes")==CMD.end())     { 
                printf("ERROR: chr sizes file is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--folder")==CMD.end()) {
                printf("ERROR: output folder  is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--marker")==CMD.end()) {
                printf("ERROR: marker    file is NOT provided. Exited.\n"); exit(1);}
            if (CMD.find("--consen")==CMD.end()) {
                printf("ERROR: consensus file is NOT provided. Exited.\n"); exit(1);}
        }
        if(CMD.find("--chromosome")!=CMD.end() && ic==argc-1) // chrmsm-region-zoom parameters
        {
            if ((CMD.find("--begin")==CMD.end() ||  CMD.find("--end")==CMD.end())){
                printf("ERROR: zoom region must be specified1. Exited.\n"); exit(1);}
            if ((CMD.find("--minfreq")==CMD.end() && CMD.find("--maxfreq")==CMD.end())){
                printf("ERROR: zoom region must be specified2. Exited.\n"); exit(1);}
            if(reg_begin>reg_end || reg_freq_min>reg_freq_max){
                printf("\nERROR: begin, end, minfreq, maxfreq are NOT proper. Exited.\n"); exit(1);}
        }
        ic ++;
    }
    
    write_log();
    write_zoom_region();
    
    return 1;
}

/*
//////  pre-check if an option is valid; if yes, put it in map //////
int precheck_opt(int argc, char* argv[], int opt_i, char* opt_name, bool arg, bool is_num)
{
    if(arg)
    {
        opt_i ++;
        if (opt_i >= argc) print_error_exit(opt_name, false);
        if(is_num)
        {
            if(!is_number(argv[opt_i])) {
                printf("%s should be a number. Exited.\n", opt_name); 
                exit(1);
            }
        }
        // record <option, arg> in map<string, string>//
        string s_key(opt_name); 
        string s_val(argv[opt_i]);
        CMD.insert(std::pair<string,string>(s_key, s_val));
    }
    else
    {
        // record <option, "true"> in map<string, string>//
        string s_key(opt_name); 
        CMD.insert(std::pair<string,string>(s_key, "true"));
    }

    return opt_i;
}

////// if there is an error on an option, exit the program //////
void print_error_exit(char* opt, bool fscanf_error)
{
    if(fscanf_error)
    {
        printf("ERROR on fscanf in %s. Exited.\n", opt);
    }
    else // cmd options
    {
        printf("ERROR: %s requires (proper) argument. Exited.\n", opt);
    }
    exit(1);
}
*/
