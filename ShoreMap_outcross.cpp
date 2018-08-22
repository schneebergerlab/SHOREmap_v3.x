/* this function visualizes frequency of individual markers and windows of markers                */
/* this function corresponds to ShoreMap_confint()  in the original implementation                */
/*  but without calculation of the confidence interval                                            */
/* Date  : 2013-March-23                                                                          */
#include                    <stddef.h>
#include                    <stdlib.h>
#include                     <stdio.h>
#include                      <math.h>
#include                      <string>
#include                    <iostream>
#include                     <sstream>
#include                    <string.h>
#include                     <ctype.h>
#include                      <time.h>
#include                         <map>
// file/dir
#include                    <dirent.h>
#include                  <sys/stat.h>
#include                 <sys/types.h>
#include                    <unistd.h>
#include           "dislin/dislin_d.h"
// sub-function
#include                   "globals.h"
#include                 "is_number.h"
#include            "cmd_initialize.h"
#include             "check_ref_err.h"
#include         "plot_chr_winboost.h"
#include        "read_allele_count2.h"
#include "filter_with_marker_parent.h"
#include               "read_marker.h"
#include     "print_filtered_marker.h"

bool visualizeFreq();

bool ShoreMap_outcross(int argc, char* argv[])
{
    /* TODO 2013-05-28 15:58: write intermediate steps into log file                              */
    if (argc < 10) {
        /* new feature filter2&3                                                 2013-06-11 19:50 */
        /*  note: --fg-N-cov and --fg-INDEL-cov can be used only when -bg-ref-filter is on        */
        printf("\nUsage: SHOREmap outcross [Option [default (if any)]].                      \n\n");
        printf("#Mandatory:                                                                    \n");
        printf("--chrsizes            STRING   File of names-and-sizes of chromosomes;   [NULL]\n");
        printf("--folder              STRING   Output folder;                            [NULL]\n");
        printf("--marker              STRING   File of markers (fg);                     [NULL]\n");
        printf("--consen              STRING   File of consensus in the pool for markers;[NULL]\n");
        printf("#Optional:                                                                     \n");
        printf("#Filter 1:                     quality criteria of alt base of markers (fg)    \n");        
        printf("--min-coverage        INT      Minimun coverage;                            [0]\n");
        printf("--max-coverage        INT      Maximum coverage;                          [INF]\n");
        printf("--marker-score        INT      Minimum quality;                             [0]\n");
        printf("--min-marker          INT      Minimum number   of markers of a window;     [2]\n");
        printf("--fg-N-cov            INT      Maximum coverage of fg 'N'-call;           [INF]\n"); 
        printf("--fg-INDEL-cov        INT      Maximum coverage of fg indel-call;         [INF]\n"); 
        printf("--marker-hit          DOUBLE   Maximum hit of a marker position (SHORE); [INF]\n");        
        ///the following are turned   off      on 2014-03-28 20:45
        ///printf("#Filter 2:                  quality criteria of ref base of markers (bg)    \n");        
        ///printf("-bg-ref-filter              Turn on filter with bg ref base info;      [off]\n");
        ///printf("--bg-ref-base-file STRING   File of bg-ref bases;                     [NULL]\n");
        ///printf("--bg-ref-cov       DOUBLE   Minimum coverage;                           [10]\n");
        ///printf("--bg-ref-cov-max   DOUBLE   Maximum coverage;                          [100]\n");
        ///printf("--bg-ref-freq      DOUBLE   Minimum concordance;                       [0.8]\n");
        ///printf("--bg-ref-score     DOUBLE   Minimum quality;                            [25]\n");
        ///printf("#Filter 3:                  bg-markers to keep or remove markers (fg)       \n");        
        ///printf("--marker-pa        STRING   File of bg markers (common kept  for fg)   [NULL]\n");
        ///printf("--marker-pb        STRING   File of bg markers (common removed from fg)[NULL]\n");
        ///printf("--pmarker-score    INT      Minimum quality   of alt base of bg-markers;  [0]\n");
        ///printf("--pmarker-min-cov  INT      Minimum coverage  of alt base of bg-markers;  [0]\n");
        ///printf("--pmarker-max-cov  INT      Maximum coverage  of alt base of bg-markers;[INF]\n");
        ///printf("--pmarker-min-freq DOUBLE   Minimum frequency of alt base of bg-markers;[0.0]\n");
        ///printf("--pmarker-ab-ratio STRING   Tune cutoffs of bg-markers;             [1,1,1,1]\n");
        printf("#Sliding-window:               analysis and visualization of AFs               \n");
        printf("--window-size         INT      Size of sliding windows;                [200000]\n");
        printf("--window-step         INT      Step of sliding windows;                 [10000]\n");
        printf("--interval-min-mean   DOUBLE   Minimum mean of AFs;                      [0.99]\n");
        printf("--interval-max-mean   DOUBLE   Maximum mean of AFs;                      [1.00]\n");
        printf("--interval-max-cvar   DOUBLE   Maximum CV of AFs;                        [0.01]\n");
        printf("-plot-boost                    Turn on  boosted-AF plotting.              [off]\n");
        printf("-plot-scale                    Turn on  scaled     plotting;              [off]\n");
        printf("-plot-win                      Turn on  plotting window-based  AFs        [off]\n");
        printf("-plot-marker-off               Turn off plotting single-marker AFs         [on]\n");
        printf("--cluster             INT      Number of clusters for ranking marker;       [5]\n");        
        printf("#Zooming:\n");           
        printf("--chromosome          STRING   ID of chromosome to zoom               [not set]\n");
        printf("--begin               INT      Begin position of zoom region          [not set]\n");
        printf("--end                 INT      End   positon  of zoom region          [not set]\n");
        printf("--minfreq             INT      Minimum AF                                 [0.0]\n");
        printf("--maxfreq             INT      Maximum AF;                                [1.0]\n");
        printf("-background2                   Mutation is in second parent;              [off]\n");
        printf("-rab                           Record point-AF, window-AF and boost values[off]\n");
        printf("-verbose                       Be talkative;                              [off]\n");
        printf("                                                                               \n");        
        printf("*bg-background;fg-foreground; AF-Allele Frequency; CV-Coefficient of Variation.\n");
        /*
        printf("*if only resequencing data of the pool is provided, we can only tune parameters\n");
        printf("    with filter 1 for filtering - mapping interval might not be accurate.      \n");
        printf("*if resequencing data of one parental line (without phenotype) is provided,    \n");
        printf("    tuning parameters with filter 2 may increase quality of markers.           \n");
        printf("    !caution: filter 2 will greatly decrease the number of markers!            \n");
        printf("*if markers are created with SHOREmap create, tunining parameters              \n");
        printf("    with filter 1 options could lead to a more accurate mapping interval.      \n");
        printf("*if there are special markers that one wants to keep or remove, use filter 3.  \n"); 
        */       
        printf("                                                                               \n");
        exit(1);
    }
    /*  defaults */
    expect              =     1.0;   // CMD["-target"]            = "";
    fchrsizes           =      "";   // CMD["--chrsizes"]         = "";
    out_folder          =      "";   // CMD["--folder"]           = "";
    fmarker             =      "";   // CMD["--marker"]           = "";
    fmarker_pa          =      "";   // CMD["--marker-pa"]        = "";
    fmarker_pb          =      "";   // CMD["--marker-pb"]        = "";
    pmarker_score       =       0;   // CMD["--pmarker-score"]    = "";
    pmarker_min_cov     =       0;   // CMD["--pmarker-min-cov"]  = "";
    pmarker_max_cov     =     INF;   // CMD["--pmarmer-max-cov"]  = "";
    pmarker_min_freq    =     0.0;   // CMD["--pmarker-min-freq"] = "";
    pmarker_ab_ratio    = "1,1,1,1"; // CMD["--pmarker-ab-ratio"] = ""; 
    /* new 2013-06-11 19:53 */
    fg_N_cov            =     INF;   // <=> marker_N_cov;     INF means no filtering with this para
    fg_INDEL_cov        =     INF;   // <=> marker_INDEL_cov; INF means no filtering with this para
    marker_hit          =(double)INF;
    bg_ref_base_file    =      "";
    bg_ref_cov          =      10;
    bg_ref_cov_max      =     100;
    bg_ref_freq         =     0.8;
    bg_ref_score        =      25;
    bg_ref_filter       =       0;
    /* end new                  */
    fconsensus          =      "";   // CMD["--consen"]            = "";
    marker_score        =       0;   // CMD["--marker-score"]      = ""; //caution: no filter if 0
    filter_min_marker   =       2;   // CMD["--min-marker"]        = "";
    filter_min_coverage =       0;   // CMD["--min-coverage"]      = ""; //caution: default no filter
    filter_max_coverage =     INF;   // CMD["--max-coverage"]      = "";

    window_size         =  200000;   // CMD["--window-size"]       = "";
    window_step         =   10000;   // CMD["--window-step"]       = "";
    interval_min_mean   =    0.99;
    interval_max_mean   =    1.00;
    interval_max_cvar   =    0.01;

    reg_chromosome      =      "";   // CMD["--chromosome"]        = "";
    reg_begin           =       0;   // CMD["--begin"]             = "";
    reg_end             =       0;   // CMD["--end"]               = "";
    reg_freq_min        =       0;   // CMD["--minfreq"]           = "";
    reg_freq_max        =       1;   // CMD["--maxfreq"]           = "";

    freferror           =      "";   // CMD["--referrors"]         = "";
    background2         =       0;   // CMD["-background2"]        = "";
    verbose             =       0;   // CMD["-verbose"]            = "";
    
    boost_max           =       0;   // CMD["-boost-max"]          = "";
    //r_max             =       0;   // CMD["-r-max"]              = "";
    plot_boost          =       0;   // related to "-plot-r"
    //plot_r            =       1;   // CMD["-plot-r"]             = "";
    plot_marker         =       1;   // CMD["-plot-marker-off"]    = "";
    plot_scale          =       0;   // CMD["-plot-scale"]         = "";
    clusterK            =       5;
                         
    runid               =       1;   // CMD["-runid"]              = "";
    /* tmp variable in this function .............................................................*/
    bool output_filtered_marker = true;
    
    /* 0.read practical values from cmd line......................................................*/
    cmd_initialize(argc, argv);  
    /* 4. write outcross log:  ...................................................................*/
    std::string fVSlog("");
    fVSlog = out_folder +  "Outcross.log";
    FILE* fplog = fopen(fVSlog.c_str(), "a+");
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
    /* new bg-ref filtering - 2013-06-11 20:08 ...................................................*/
    /* 1.5.filter foreground SNPs (candidate markers) if option is set............................*/
    if(bg_ref_filter == 1)
    {
        std::string fg_snp_file = fmarker;
        std::string fg_snp_file_filtered = "";           // markers satisfying filtering constraints
        map<std::string, MARK6> fg_snp_map_filtered;
        check_ref_err((char*)fg_snp_file.c_str(),
                      (char*)fconsensus.c_str(), 
                      (char*)bg_ref_base_file.c_str(), 
                      bg_ref_cov, 
                      bg_ref_cov_max,
                      bg_ref_freq,
                      fg_N_cov, 
                      fg_INDEL_cov, 
                      bg_ref_score,
                      &fg_snp_file_filtered,
                      &fg_snp_map_filtered);
        fmarker.clear();
        fmarker = fg_snp_file_filtered;                 // set new set of SNPs as fmarker...........
        fprintf(fplog,"\nFG markers in %s (filtered according to BG-ref-base quality, concordance):",
                fg_snp_file.c_str());
        fprintf(fplog, "\n\t\t%s\n", fg_snp_file_filtered.c_str());
    }
    /* normal filtering */
    /* 1.read bg-ref-filtered markers: ALLELE1.key == ALLELE2.key == chr+".#."+pos................*/
    unsigned long num_m_given = 0;
    unsigned long num_m_inuse = 0;
    if(!read_marker((char*)fmarker.c_str(), &num_m_given, &num_m_inuse))
    {
       printf("ERROR: no marker info recorded. Exited. \n"); exit(1);
    }
    /* 1.5 filter markers with SNPs from background A and/or B (if sequenced).....................*/
    /* idea: assume bg A with the same phenotype as the pool while bg B without.....................
             After calling SNPs for A/B (and the pool) from resequencing against a reference seq...: 
             case 1: keep   pool-SNPs present in A.................................................. 
             case 2: remove pool-SNPs present in B..................................................
             NOTE: Quality of SNPs of A/B can be checked before using them..........................
    */
    if(fmarker_pa.length() > 0)
    {
        /* SNPs of the bg with the same phenotype.................................................*/
        if(verbose) 
        {    
            printf("Filtering markers in the pool by using variations of a background; ");
            printf("this background has the same phenotype as the mutant.\n");
            printf("Number of pool-markers we now have: %ld\n", ALLELE2.size());
        }
        if(!filter_with_marker_parent((char*)fmarker_pa.c_str(), true))
        {
            printf("No common variations found between background A and the pool. \n");
            printf("Pls find a proper SNP file if you want to filter with --marker-pa;\n");
            printf("Otherwise, remove option --marker-pa from the cmd line. Exited. \n");
            exit(1);
        }
        if(verbose) 
        {  
            printf("Number of pool-markers after filtering with this parent: %ld.\n", 
                   ALLELE2.size());
        }
    }
    if(fmarker_pb.length() > 0)
    {
        /* SNPs of the bg without the phenotype...................................................*/
        if(verbose) 
        {    
            printf("Filtering markers in the pool by using variations of background; ");
            printf("this background does not have the interested phenotype of the mutant.\n");
            printf("Number of pool-markers we now have: %ld\n", ALLELE2.size());
        }
        if(!filter_with_marker_parent((char*)fmarker_pb.c_str(), false))
        {
            printf("All markers of the pool have been removed according to background B. \n");
            printf("Pls find a proper SNPs file if you want to filter with --marker-pb;\n");
            printf("Otherwise, remove option --marker-pb from the cmd line. Exited. \n");
            exit(1);
        }
        if(verbose) 
        {  
            printf("Number of pool-markers different from this background: %ld.\n", 
                   ALLELE2.size());
        }
    }
    
    /* 2. read consensus base and allele1 allele2 error counts of mutant-markers..................*/
    if (!read_allele_counts2((char*)fconsensus.c_str()))
    {
        printf("ERROR: no consensus info recorded. Exited. \n"); 
        exit(1);
    }
    
    /* 3. visualize frequency of markers of mutant-markers........................................*/
    /*    note: CHR2POS2_ale1_ale2_err_COUNT can be filtered with marker freq, cov, etc...........*/
    if(!visualizeFreq())
    {
        printf("ERROR: visualization failed. Exited. \n"); 
        exit(1);
    }
    
    /* 3.5output markers to visualize.............................................................*/
    if(output_filtered_marker)
    {
        map<std::string, map<unsigned long, TRIPLE> > tCNT;
        unsigned long num_marker_consen = 0;
        map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr;
        map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr_end;
        chr_itr     = CHR2POS2_ale1_ale2_err_COUNT.begin();
        chr_itr_end = CHR2POS2_ale1_ale2_err_COUNT.end();
        while(chr_itr != chr_itr_end)
        {
            if(CMD.find("--chromosome")!=CMD.end())
            {
                if(reg_chromosome==(*chr_itr).first)
                {
                    tCNT[(*chr_itr).first].insert((*chr_itr).second.begin(), (*chr_itr).second.end());
                    num_marker_consen += (*chr_itr).second.size();
                    break;
                }
                else
                {
                    chr_itr ++;
                    continue;
                }
            }
            else
            {
                num_marker_consen += (*chr_itr).second.size();
                chr_itr ++;
            }
        }
        if(verbose) printf("Recording consen info of %ld markers after visualization. \n",
                           num_marker_consen);
        if(CMD.find("--chromosome") == CMD.end())
        {
            tCNT.insert(CHR2POS2_ale1_ale2_err_COUNT.begin(), CHR2POS2_ale1_ale2_err_COUNT.end());
        }
        if(print_filtered_marker(tCNT, ALLELE2, ALLELE1, QUALITY1))
        {
            fprintf(fplog, "Markers are further filtered ");
            fprintf(fplog, "with user-defined/default parameters, and recorded. \n");
        }
    }
    fprintf(fplog, "\n\nOutputs achieved.\n");
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    fprintf(fplog, "\nOutcross analysis and visualization are successfully finished on %s\n", 
            asctime(tinfo));
    fclose(fplog);
    
    return true;
}

bool visualizeFreq()
{
    double        foreground_frequency      = expect;
    map<std::string, unsigned long> chrsize = CHR2SIZE;
    unsigned long filterOutliers    = outlier_window_size;
    double        filterPValue      = outlier_pvalue;
    unsigned long minMarker         = filter_min_marker;
    unsigned long minCoverage       = filter_min_coverage;
    unsigned long maxCoverage       = filter_max_coverage;
    //double      boostMax          = boost_max;
    //bool        plotBoost         = plot_boost;
    //bool        plotMarker        = plot_marker;
    unsigned long runID             = runid;    
    
    if(false && verbose) // TO REMOVE
    {
        printf("Chromosome info for SHOREmap visualize: \n");
        map<std::string, unsigned long>::iterator iter;
        for (iter = chrsize.begin(); iter != chrsize.end(); iter++) 
        {
            printf("\t\t\t\t\t%s, %ld\n", (*iter).first.c_str(), (*iter).second);
        }
    }
    /* set up zoom interval if it exists: no need to set, just use them directly ----- 2013-03-28 */
    if(CMD.find("--chromosome")!=CMD.end())
    {
        std::string   z_chr = reg_chromosome;
        unsigned long z_beg = reg_begin;
        unsigned long z_end = reg_end;
        double        z_min = reg_freq_min;
        double        z_max = reg_freq_max;
    }
    
    /* begin: visualize the chromosomes one by one................................................*/
    
    /* set up pdf file */
    char chrNameBuf[512];
    if(plot_boost)
    sprintf(chrNameBuf, "%sOC_AF_visualization_%ld_boost_%s.pdf\0", out_folder.c_str(), runid, reg_chromosome.c_str());
    else
    sprintf(chrNameBuf, "%sOC_AF_visualization_%ld_%s.pdf\0",       out_folder.c_str(), runid, reg_chromosome.c_str());
    metafl("PDF");                         // level 0.....- defines the metafile format.............
    setfil(chrNameBuf);                    // level 0.....- sets alternate filename.................
    filmod("VERSION");                     // level 0/1/2/3 do not shorten existing filename........
    scrmod("revers");                      // level 0.....- swaps back and foreground colours.......
    setpag("USEL");                        // level 0.....- selects a predefined page format........
                                           //        ...... landscape: 11180*8640 points............
    /* initialization                      */
    disini();                              // level 0.....- initialize DISLIN with default PARAs and
                                           //        ...... a plotfile. The level is set to 1. .....
                                           //        ...... DISINI must be called before any other.. 
                                           //        ...... DISLIN routine except for ..............
    unsigned long plot_ith_chr = 0;
    map<std::string, unsigned long>::iterator chr_itr;
    for(chr_itr = chrsize.begin(); chr_itr != chrsize.end(); chr_itr ++)
    {
        std::string ci_chromosome = (*chr_itr).first;
        if(CMD.find("--chromosome")!=CMD.end() && ci_chromosome!=reg_chromosome) 
        {
            /* check 1: skip chromosomes until finding the specified one..........................*/
            continue;
        }
        if(CHR2SIZE.find(ci_chromosome) == CHR2SIZE.end())
        {
            /* check 2: if a chromosome is in the interested list.................................*/
            printf("Chromosome %s not found in the given list, skipped.\n\n", ci_chromosome.c_str());
            continue;
        }
        /* print info on screen */
        printf("Start visualizing chromosome %s", ci_chromosome.c_str());
        if(CMD.find("--chromosome") != CMD.end())
        {
            printf(" (as specified)...\n");
            /* parameters: reg_begin; reg_end; reg_freq_min; reg_freq_max;........................*/
            /* will be used for filtering if a chr is specified...................................*/
        }
        else
        {
            printf("...\n");
        }
        /* get counts of allele 1, allele2, error for all markers of a single chromosome ci.......*/
        map<unsigned long, TRIPLE> internalData = CHR2POS2_ale1_ale2_err_COUNT[(*chr_itr).first];
        if(internalData.size() == 0) 
        {
            printf("\tNo markers found for %s. Visualization skipped. \n\n", 
                  (char*)ci_chromosome.c_str());
            continue;
        }
        
        /* begin: print info for checking if necessary............................................*/
        if(false && verbose)
        {
            printf("\tnumber of markers recorded: %ld\n", internalData.size());
            map<unsigned long, TRIPLE>::iterator it2;
            for (it2 = internalData.begin(); it2 != internalData.end(); it2++)
            {
                printf("\t\t%s.#.%ld: \t", ci_chromosome.c_str(), (*it2).first);
                printf("%ld\t",   (*it2).second.Ci[0]);
                printf("%ld\t",   (*it2).second.Ci[1]);
                printf("%ld\t\n", (*it2).second.Ci[2]);                           
            }
        }
        /* end: print info for checking if necessary..............................................*/
    
        /* begin: 1st filtering: coverage (if specificed, with region range, allele frequency)....*/
        unsigned long num_record_b4_1st_filter = internalData.size();
        if(verbose) 
        {
            printf("\tnumber of markers before filtering with coverage: %ld.\n", 
                   num_record_b4_1st_filter);
        }
        map<unsigned long, TRIPLE>::iterator pos_itr;
        pos_itr = internalData.begin();
        while(pos_itr != internalData.end())
        {
            unsigned long ipos = (*pos_itr).first;
            TRIPLE rcod        = (*pos_itr).second;
            //unsigned long cvrg = rcod.Ci[0] + rcod.Ci[1] + rcod.Ci[2];  // caution: error  counted
            unsigned long cvrg = rcod.Ci[0] + rcod.Ci[1];                 // caution: error not used
            bool is_in_REG     = true;
            bool is_gd_AFQ     = true;
            bool is_gd_COV     = true;  
             /* check if in interested region.....................................................*/
            if(CMD.find("--chromosome") != CMD.end())
            {
                if(ipos>reg_end || ipos<reg_begin)           
                {
                    is_in_REG = false;
                }
                else
                {
                /* check if satisfactory allele frequency   */
                    //double iAF = (double)rcod.Ci[1]/cvrg;
                    //if(iAF>reg_freq_max || iAF<reg_freq_min) 
                    //{
                    //    is_gd_AFQ = false;
                    //}
                }
            }
            /* check if satisfactory allele frequency.............................................*/
            double iAF = (double)rcod.Ci[1]/cvrg;
            if(iAF>reg_freq_max || iAF<reg_freq_min) 
            {
                is_gd_AFQ = false;
            }
            /* check if satisfactory coverage.....................................................*/
            bool icov_flag = true;
            std::stringstream ss;
            ss << ipos;
            std::string fndkey = (ci_chromosome + ".#." + ss.str());
            
            if(is_in_REG && is_gd_AFQ)
            {
                if(QUALITY1[fndkey].find("FLAG4parentB") != std::string::npos)
                {
                    /* 2.13-09-20> values of bg-a & b were swapped when reading consensus info... */
                    /* so we have to check the other value which is for the current marker....... */
                    //if( (rcod.Ci[0]>maxCoverage) || (rcod.Ci[0]<minCoverage)) // caution here
                    //{
                    //    is_gd_COV = false;
                    //}
                    /* normal coverage as swap is done when reading consen     - 2014-04-23 18:41 */
                    if( (rcod.Ci[1]>maxCoverage) || (rcod.Ci[1]<minCoverage)) // caution here
                    {
                        is_gd_COV = false;
                    }
                }
                else
                {
                    /* normal case */
                    if( (rcod.Ci[1]>maxCoverage) || (rcod.Ci[1]<minCoverage))
                    {
                        is_gd_COV = false;
                    }            
                }
            }
            
            if(!is_in_REG || !is_gd_AFQ || !is_gd_COV)
            {
                internalData.erase(pos_itr++);
            }
            else pos_itr ++;
        }
        if(verbose) 
        {
            printf("\t..................after.........................: %ld.\n", 
                   internalData.size());
        }
        CHR2POS2_ale1_ale2_err_COUNT[(*chr_itr).first].clear();
        if(internalData.size() == 0) 
        {
            printf("\tNo markers left for %s after filtering with coverage. ",
                   (char*)ci_chromosome.c_str());     
            printf("Visualization skipped. \n\n");
            continue;
        }
        /* update consen records if there is records in internalData..............................*/ 
        CHR2POS2_ale1_ale2_err_COUNT[(*chr_itr).first].insert(internalData.begin(), internalData.end());
        /* end: 1st filtering with coverage */
        map<unsigned long, TRIPLE> filtered;
        filtered.clear();
        /* second-time filetering of outliers has not been used   -  2013-03-23 */
        // filterSampling(&internalData, filterOutliers, filterPValue, &filtered);
        
        /* begin: visualizing chromosome   */
        plot_ith_chr ++;                   // the ith chromosome to plot............................
        if(plot_ith_chr > 1) newpag();
        plot_chr_winboost(ci_chromosome, plot_ith_chr, internalData, filtered);
        printf("Visualization of chromosome %s done. \n\n", ci_chromosome.c_str());
        /* end:   visualizing chromosome   */
    }
    /* end:  visualizing the chromosomes one by one...............................................*/   
    /* close pdf file */ 
    disfin();                              // level 1/2/3 - terminates DISLIN and prints a message on
                                           //             ..the screen. The level is set back to 0.
                                           //             ..not required for qplsca(...)
    
    return true;
}
