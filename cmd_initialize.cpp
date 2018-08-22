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

/* initialize outcross parameters from cmd line */

int  cmd_initialize(int argc, char* argv[])
{
    /* 
      Read parameters from cmd options; option-values are checked and updated in ...................
      map<string, string> CMD with format: CMD["option"] = "arg" (if exist).........................
    */
    int ic = 2;                               // check verbose first................................
    while (ic < argc)
    {
        if(!strcmp(argv[ic],"-verbose"))      // option 27 to variable: verbose.....................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            verbose         = 1;
            printf("Be talkative during process. \n");
        }
        ic ++;
    }
    ic = 2;    // option ic=0: SHOREmap; ic=1: outcross/backcross/extract/annotate..................
    while (ic < argc) 
    {
        strcatCMD += " ";
        strcatCMD += (string)argv[ic];
        if(!strcmp(argv[ic],"--chrsizes"))    // option 1 to variable: fchrsizes....................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fchrsizes        += argv[ic];     // set value for option...............................
            FILE* fp_chrsizes = fopen((char*)fchrsizes.c_str(), "r");
            if(!fp_chrsizes)  { print_error_exit(argv[ic-1], false); }
            fclose(fp_chrsizes);              // check done.........................................
            if (verbose) printf("Chrs sizes read from file:\t\t%s\n", fchrsizes.c_str());
            if(!read_chromosomes((char*)fchrsizes.c_str()))  // read contents.......................
            {
                printf("ERROR: invalid content in file %s.\n", fchrsizes.c_str());
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--folder")) // option 2 to variable: out_folder...................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            out_folder         += argv[ic];
            if(out_folder[out_folder.length()-1] != '/') out_folder += "/";
            DIR* dir_out_folder = opendir((char*)out_folder.c_str());
            if(dir_out_folder  == NULL)
            {
                if(!mkdir((char*)out_folder.c_str(), S_IRWXU|S_IRWXG|S_IRWXO))
                {
                 /* if !mkdir() is TRUE: a new directory has to be created..........................
                    if out_folder is "outfolder/" instead of "/your/output/path/outfolder".......... 
                    an "outfolder/" will be created under the current working directory...........*/
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
        else if(!strcmp(argv[ic],"--marker"))     // option 3 to variable: fmarker..................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker        += argv[ic];
            FILE* fp_marker = fopen((char*)fmarker.c_str(), "r"); 
            if(fp_marker == NULL)
            {
                printf("snp file \'%s\' does NOT exist. Exited.\n", fmarker.c_str());
                exit(1);
            }
            if (verbose) printf("File of snps provided:\t\t\t%s\n", fmarker.c_str());
           // read_marker();
            fclose(fp_marker);
        }
        else if(!strcmp(argv[ic],"--marker-pa"))  // option 3.1 to variable: fmarker_pa.............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker_pa        += argv[ic];
            FILE* fp_marker_pa = fopen((char*)fmarker_pa.c_str(), "r"); 
            if(fp_marker_pa == NULL)
            {
                printf("snp-file of parent a \'%s\' does NOT exist. Exited.\n", fmarker_pa.c_str());
                exit(1);
            }
            if (verbose) printf("File of snps provided as parent a:\t%s\n", fmarker_pa.c_str());
            fclose(fp_marker_pa);
        }
        else if(!strcmp(argv[ic],"--marker-pb"))  // option 3.2 to variable: fmarker_pb.............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker_pb        += argv[ic];
            FILE* fp_marker_pb = fopen((char*)fmarker_pb.c_str(), "r"); 
            if(fp_marker_pb == NULL)
            {
                printf("snp-file of parent b \'%s\' does NOT exist. Exited.\n", fmarker_pb.c_str());
                exit(1);
            }
            if (verbose) printf("File of snps provided as parent b:\t%s\n", fmarker_pb.c_str());
            fclose(fp_marker_pb);
        }
        else if(!strcmp(argv[ic],"--consen"))    // option 4 to variable: fconsensus................
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
        else if(!strcmp(argv[ic],"--interval-min-mean")) // option 5 to variable: interval_min_mean.
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_min_mean       = atof(argv[ic]);
            if (verbose) 
            {
               printf("interval_min_mean set as:\t\t%.4f\n", 
                      interval_min_mean);
            }
        }
        else if(!strcmp(argv[ic],"--interval-max-mean"))// option 5.5 to variable: interval_max_mean
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_max_mean       = atof(argv[ic]);
            if (verbose) 
            {
               printf("interval_max_mean set as:\t\t%.4f\n", 
                      interval_max_mean);
            }
        }
        else if(!strcmp(argv[ic],"--interval-max-cvar"))// option 6 to variable: interval_max_cvar..
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_max_cvar    = atof(argv[ic]);
            if (verbose)
            {
                printf("interval_max_cvar set as:\t\t%.4f.\n", 
                        interval_max_cvar);
            }
        }
        else if(!strcmp(argv[ic],"--peak-window-size")) // option 7 to variable: peak_window_size...
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            peak_window_size          = atol(argv[ic]);
            if(peak_window_size<1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("peak_window_size set as:\t\t%ld\n", peak_window_size);
        }
        else if(!strcmp(argv[ic],"--peak-window-step")) // option 8 to variable: peak_window_step...
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            peak_window_step          = atol(argv[ic]);
            if(peak_window_step<1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("peak_window_step set as:\t\t%ld\n", peak_window_step);
        }
        else if(!strcmp(argv[ic],"--mis-phenotyped"))  // option 9 to variable: misphenotyped.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            misphenotyped           = atof(argv[ic]);
            if(misphenotyped<0.0 || misphenotyped>1.0) 
                { printf("ERROR: arg of %s must be [0.0, 1.0].\n", argv[ic-1]); exit(1);}
            if (verbose) printf("misphenotyped set as:\t\t\t%.2f\n", misphenotyped);
        }
        else if(!strcmp(argv[ic],"--window-size"))     // option 10 to variable: window_size........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_size           = atol(argv[ic]);
            if(window_size<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("window_size set as:\t\t\t%ld\n", window_size);
        }
        else if(!strcmp(argv[ic],"--window-step"))     // option 11 to variable: window_step........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_step           = atol(argv[ic]);
            if(window_step<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("window_step set as:\t\t\t%ld\n", window_step);
        }
        else if(!strcmp(argv[ic],"-plot-marker-off"))  // option 12 to variable: plot_marker........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_marker         = 0;
            if (verbose) printf("plot_marker set as:\t\t\t%d (off marker-af-plot)\n", plot_marker);
        }
        else if(!strcmp(argv[ic],"-plot-boost"))       // option 13 to variable: plot_r.............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_r         = 0;
            plot_boost     = 1; // caution
            if (verbose) printf("plot_boost set as:\t\t\t%d\n", plot_boost);
        }
        else if(!strcmp(argv[ic],"-plot-win"))         // option 22.5 to variable: plot_window......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_window     = 1; // caution
            if (verbose) printf("plot_window set as:\t\t\t%d\n", plot_window);
        }
        else if(!strcmp(argv[ic],"-plot-scale"))       // option 13.5 to variable: plot_scale.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_scale     = 1; // caution
            if (verbose) printf("plot_scale set as:\t\t\ttrue\n", plot_scale);
        }
        else if(!strcmp(argv[ic],"--min-marker"))      // option 14 to variable: filter_min_marker..
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_marker   = atol(argv[ic]);
            if(filter_min_marker<=1)
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("filter_min_marker set as:\t\t%ld\n", filter_min_marker);
        }
        else if(!strcmp(argv[ic],"--min-coverage"))    // option 15 to variable: filter_min_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_min_coverage set as:\t\t%ld\n", filter_min_coverage);
        }
        else if(!strcmp(argv[ic],"--max-coverage"))    // option 16 to variable: filter_max_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_max_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_max_coverage set as:\t\t%ld\n", filter_max_coverage);
        }
        else if(!strcmp(argv[ic],"--marker-hit"))      // option 16.4 to variable: marker_hit.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            marker_hit           = atof(argv[ic]);
            if(marker_hit<0) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("marker_hit:\t\t\t\t%.4f\n", marker_hit);
        }
        else if(!strcmp(argv[ic],"--marker-score"))    // option 16.5 to variable: marker_score.....
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            marker_score           = atol(argv[ic]);
            if(marker_score<=1)
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("marker_score:\t\t\t\t%ld\n", marker_score);
        }
        else if(!strcmp(argv[ic],"--pmarker-score"))   // option 16.6 to variable: pmarker_score....
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_score           = atol(argv[ic]);
            if(pmarker_score < 1)
                { printf("ERROR: arg of %s must be >= 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("pmarker_score:\t\t\t%ld\n", pmarker_score);
        }
        else if(!strcmp(argv[ic],"--pmarker-min-cov")) // option 16.7 to variable: pmarker_min_cov..
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_min_cov   = atol(argv[ic]);
            if (verbose) printf("pmarker_min_cov set as:\t\t\t%ld\n", pmarker_min_cov);
        }
        else if(!strcmp(argv[ic],"--pmarker-max-cov")) // option 16.8 to variable: pmarker_max_cov..
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_max_cov   = atol(argv[ic]);
            if (verbose) printf("pmarker_max_cov set as:\t\t\t%ld\n", pmarker_max_cov);
        }
        else if(!strcmp(argv[ic],"--pmarker-min-freq"))// option 16.9 to variable: pmarker_min_freq.
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pmarker_min_freq     = atof(argv[ic]);
            if(pmarker_min_freq<0.0 || pmarker_min_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            if (verbose) printf("pmarker_min_freq set as:\t\t%.2f.\n", pmarker_min_freq);
        }
        else if(!strcmp(argv[ic],"--pmarker-ab-ratio"))// option 16.95 to variable: pmarker_ab_ratio
        {
            /* suppose cutoffs are for parent a */
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            pmarker_ab_ratio = argv[ic];
            if (verbose) 
            {
                printf("quality-score, min-coverage, max-coverage, min-frequency = ");
                printf("Cutoff ratio of parents a,b:\t\t%s\n", pmarker_ab_ratio.c_str());
            }
        }
        else if(!strcmp(argv[ic],"-outlier"))          // option 17 to variable: outlier_window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            outlier_window_size = 0;
            if (verbose) printf("outlier removal is on. \n");
        }
        else if(!strcmp(argv[ic],"--outlier-window-size"))// option 18 to variable: not used!!!!!!!!
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            outlier_window_size          = atol(argv[ic]);
            if (verbose) printf("outlier_window_size set as:\t\t%ld\n", outlier_window_size);
        }
        else if(!strcmp(argv[ic],"--outlier-pvalue"))   // option 19 to variable: outlier_pvalue....
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            outlier_pvalue          = atof(argv[ic]);
            if(outlier_pvalue<0.0 || outlier_pvalue>1.0)
                { printf("ERROR: arg of %s must be set in [0.0, 1.0]. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("outlier_pvalue set as:\t\t\t%.2f.\n", outlier_pvalue);
        }
        else if(!strcmp(argv[ic], "--chromosome"))      // option 20 to variable: reg_chromosome....
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
        else if(!strcmp(argv[ic],"--begin"))            // option 21 to variable: reg_begin.........
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
        else if(!strcmp(argv[ic],"--end"))              // option 22 to variable: reg_end...........
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
        else if(!strcmp(argv[ic],"--minfreq"))          // option 23 to variable: reg_freq_min......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_freq_min     = atof(argv[ic]);
            if(reg_freq_min<0.0 || reg_freq_min>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            if (verbose) printf("reg_freq_min set as:\t\t\t%.2f.\n", reg_freq_min);
        }
        else if(!strcmp(argv[ic],"--maxfreq"))          // option 24 to variable: reg_freq_max......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_freq_max     = atof(argv[ic]);
            if(reg_freq_max<0.0 || reg_freq_max>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("reg_freq_max set as:\t\t\t%.2f.\n", reg_freq_max);
        }
        else if(!strcmp(argv[ic],"--referrors"))        // option 25 to variable: freferror.........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            freferror         += argv[ic];
            FILE* fp_referror  = fopen((char*)freferror.c_str(), "r"); 
            if(fp_referror == NULL)
            {printf("Ref-error file \'%s\' does NOT exist. Exited.\n", freferror.c_str()); exit(1);}
            fclose(fp_referror);
        }
        else if(!strcmp(argv[ic],"-background2"))       // option 26 to variable: background2.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            background2 = 1;
            if (verbose) 
            {
                printf("background2 set as:\t\t\t1 ");
                printf("(Calculation of AFs of first-column alleles has been told) \n");
            }
        }
        else if(!strcmp(argv[ic],"-runid"))             // option 27 to variable: runid.............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            runid = atoi(argv[ic]);
            if (verbose) printf("runid set as:\t\t\t%ld.\n", runid);
        }
        else if(!strcmp(argv[ic],"-boost-max"))         // option 28 to variable: boost_max.........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            boost_max         = atof(argv[ic]);
            if (verbose) printf("boost_max has set as:\t\t\t%.2f.\n", boost_max);
        }
        else if(!strcmp(argv[ic],"-r-max"))             // option 29 to variable: r_max.............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            r_max         = atof(argv[ic]);
            if (verbose) printf("r_max has set as:\t\t\t%ld.\n", r_max);
        }
        else if(!strcmp(argv[ic],"-target"))            // option 30 to variable: expect............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            expect         = atof(argv[ic]);
            if (verbose) printf("target/expect has set as:\t\t%.2f.\n", expect);
        }
        else if(!strcmp(argv[ic],"-bg-ref-filter"))     // option 31 to variable: bg_ref_filter.....
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            bg_ref_filter = 1;
            if (verbose) printf("bg_ref_filter set as:\t\t%ld\n", bg_ref_filter);
        }
        else if(!strcmp(argv[ic],"--fg-N-cov"))         // option 32 to variable: fg_N_cov..........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            fg_N_cov   = atol(argv[ic]);
            if (verbose) printf("fg_N_cov set as:\t\t\t%ld\n", fg_N_cov);
        }
        else if(!strcmp(argv[ic],"--fg-INDEL-cov"))     // option 33 to variable: fg_INDEL_cov......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            fg_INDEL_cov   = atol(argv[ic]);
            if (verbose) printf("fg_INDEL_cov set as:\t\t\t%ld\n", fg_INDEL_cov);
        }
        else if(!strcmp(argv[ic],"--bg-ref-base-file")) // option 34 to variable: bg_ref_base_file..
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
        else if(!strcmp(argv[ic],"--bg-ref-cov"))       // option 35 to variable: bg_ref_cov........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov set as:\t\t\t%ld\n", bg_ref_cov);
        }
        else if(!strcmp(argv[ic],"--bg-ref-cov-max"))   // option 35.5 to variable: bg_ref_cov_max..
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov_max   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov_max set as:\t\t\t%ld\n", bg_ref_cov_max);
        }
        else if(!strcmp(argv[ic],"--bg-ref-freq"))      // option 36 to variable: bg_ref_freq.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_freq     = atof(argv[ic]);
            if(bg_ref_freq<0.0 || bg_ref_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("bg_ref_freq:\t\t\t\t%.2f.\n", bg_ref_freq);
        }
        else if(!strcmp(argv[ic],"--bg-ref-score"))     // option 37 to variable: bg_ref_score......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_score           = atol(argv[ic]);
            if(bg_ref_score<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("bg_ref_score:\t\t\t\t%ld\n", bg_ref_score);
        }
        else if(!strcmp(argv[ic],"-rab"))               // option 38 to variable: plot_record.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_record     = 1;
            if (verbose) printf("plot_record set as:\t\t\t%ld\n", plot_record);
        }
        else if(!strcmp(argv[ic],"--cluster"))          // option 39 to variable: clusterK..........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            clusterK   = atol(argv[ic]);
            if(clusterK<1 || clusterK>40) 
                { printf("ERROR: arg of %s must be in [1, 40]. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("clusterK set as:\t\t\t%ld\n", clusterK);
        }
        else // other options not necessary.........................................................
        { 
             if(!strcmp(argv[ic],"-verbose"));          // this has been set........................
             else printf("Warning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        
        if(ic == argc-1) // check mandatory options.................................................
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
        if(CMD.find("--interval-min-mean")!=CMD.end() || CMD.find("--interval-max-mean")!=CMD.end())
        {
            if(interval_min_mean > interval_max_mean){
                printf("ERROR: interval-min-mean > interval-max-mean. Exited.\n"); exit(1);}
        }
        if(CMD.find("--chromosome")!=CMD.end() && ic==argc-1) // chrmsm-region-zoom parameters......
        {
            if ((CMD.find("--begin")==CMD.end() ||  CMD.find("--end")==CMD.end())){
                printf("ERROR: zoom region must be specified 1. Exited.\n"); exit(1);}
            if ((CMD.find("--minfreq")==CMD.end() && CMD.find("--maxfreq")==CMD.end())){
                printf("ERROR: zoom region must be specified 2. Exited.\n"); exit(1);}
            if(reg_begin>reg_end || reg_freq_min>reg_freq_max){
                printf("ERROR: begin, end, minfreq, maxfreq are NOT proper. Exited.\n"); exit(1);}
        }
        ic ++;
    }
    
    cluster_avg_coverage = (filter_min_coverage + filter_max_coverage)/2;
    cluster_dim_coverage = (filter_max_coverage - filter_min_coverage)/2;
    if(cluster_dim_coverage == 0) cluster_dim_coverage = 1; // caution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(clusterK >= 2)
    if(filter_max_coverage==INF && filter_min_coverage==0)
    {
        printf("WARNING: there is not any constraint on minimum and maximum coverage of mut base.\n");
        printf("	 check option --max-coverage for a better clustering of markers. \n");
        printf("	 average coverage for clusterting will be assumed as 50, ");
        printf("	 e.g., coverage [25, 50, 75] will receive scaling as [0.5, 1, 0.5].\n");
        printf("	 if 50 is far different from your data, please provide it. \n");
        cluster_avg_coverage = 50;
        cluster_dim_coverage = 25;
    }
    else if(filter_max_coverage==INF)
    {
        printf("WARNING: there is only constraint on minimum coverage of mut base.\n");
        printf("	 check option --max-coverage for a better clustering of markers. \n");
        printf("	 average coverage for clusterting will be assumed as 2*min-coverage, ");
        printf("e.g., coverage [%ld, %ld, %ld] will receive scaling as [0.5, 1, 0.5]. \n",
               filter_min_coverage, 2*filter_min_coverage, 3*filter_min_coverage);
        printf("	 if this average is far different from your data, please provide it. \n");
        cluster_avg_coverage = filter_min_coverage*2;
        cluster_dim_coverage = filter_min_coverage;
    }
    
    return 1;
}
