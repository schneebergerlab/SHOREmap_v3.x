#include             <stddef.h>
#include             <stdlib.h>
#include              <stdio.h>
#include               <math.h>
#include              <ctype.h>
#include             <string.h>
#include               <time.h>
#include                  <map>
#include               <vector>
#include               <string>
// file, dir
#include             <dirent.h>
#include           <sys/stat.h>
#include          <sys/types.h>
#include             <unistd.h>
// visualization
#include    "dislin/dislin_d.h"
#include            "globals.h"
#include          "is_number.h"
#include       "precheck_opt.h"
#include   "print_error_exit.h"
#include       "split_string.h"
#include   "read_chromosomes.h"
#include          "filter_fg.h"
#include      "check_ref_err.h"
#include        "read_marker.h"
//#include      "plot_chr_AF.h"
#include  "plot_chr_winboost.h"
#include "read_allele_count2.h"

using namespace std;
// sub  function
/* initialize parameters according to the cmd line input    */
int cmd_backcross(int argc, char* argv[]);
/* clean foreground markers according to background markers */
int filter_fg_with_bg(char* fileBG);
/* visualize allele frequency                               */
bool visualizeFreq_bc();

// main function
int init_backcross(int argc, char* argv[])
{
    if(argc < 3)
    {
        /* Backcross selects SNP markers with foreground/pool and background/parental info. Paras:*/
        /* "marker/fg": info on interested SNPs; "bg": info on SNPs for filtering SNPs of pool,   */
        /*  where "bg" info can be: "--bg-*", or "--bg-ref-*".                                    */
        /*  note:  chromosomes not in chrsizes file also excluded (filter_fg_with_bg(...)).       */
        /*  note2: --fg-N-cov and --fg-INDEL-cov can be used only when -bg-ref-filter is on       */
        printf("\nUsage: SHOREmap backcross [option [default (if any)]]                     \n\n");
        printf("#Mandatory:                                                                   \n");
        printf("--chrsizes            STRING   File of names and sizes of chromosomes;  [NULL]\n");
        printf("--folder              STRING   Output folder;                           [NULL]\n");
        printf("--marker              STRING   File of SNPs (F2/ref) as candidate markers;[NULL]\n");
        printf("#Optional:                                                                    \n");
        printf("#Filter 1:                     quality criteria of alt base of markers (fg)   \n");        
        printf("--marker-score        INT      Minimum quality;                           [25]\n"); 
        printf("--marker-freq         INT      Minimum concordances;                     [0.0]\n");
        printf("--min-coverage        INT      Minimun coverage;                           [1]\n");
        printf("--max-coverage        INT      Maximum coverage;                         [INF]\n");
        printf("--marker-hit          INT      Maximum hit of a marker position (SHORE); [INF]\n");
        printf("--fg-N-cov            INT      Maximum coverage of fg 'N'-call (SHORE);  [INF]\n"); 
        printf("--fg-INDEL-cov        INT      Maximum coverage of fg indel-call (SHORE);[INF]\n");
        printf("#Filter 2:                     use bg-SNPs to filter fg-markers               \n");     
        printf("--bg                  STRING   File of bg mutations;                    [NULL]\n");
        printf("--bg-score            INT      Minumum score;                              [1]\n");
        printf("--bg-freq             INT      Minumum concordance;                      [0.0]\n");
        printf("--bg-cov              INT      Minimum read support;                       [1]\n");
        printf("#Filter 3:                     quality criteria of ref base of markers (fg)   \n");        
        printf("-bg-ref-filter                 Turn on filter with bg ref base info;     [off]\n");
        printf("--consen              STRING   File of fg-consensus; necessary if plotting;[NULL]\n");
        printf("--bg-ref-base-file    STRING   File of ref bases   of bg-SNPs;          [NULL]\n");
        printf("--bg-ref-cov          DOUBLE   Minimum coverage    of bg ref base;        [10]\n");
        printf("--bg-ref-cov-max      DOUBLE   Maximum coverage    of bg ref base;       [500]\n");
        printf("--bg-ref-freq         DOUBLE   Minimum concordance of bg ref base;       [0.8]\n");
        printf("--bg-ref-score        DOUBLE   Minimum quality     of bg ref base;        [25]\n");
      //printf("#Sliding-window:               Analysis and visualization                     \n");
        printf("-plot-bc                       Turn on plotting AF;                      [off]\n");
        printf("-plot-win                      Turn on plotting window-based AFs         [off]\n");
        printf("-plot-scale                    Turn on  scaled     plotting;              [off]\n");        
        printf("--cluster             INT      Number of clusters for ranking marker;     [10]\n");
      //printf("-plot-boost                    Turn on plotting boosted AF;              [off]\n");
        printf("--window-size         INT      Size of sliding windows;               [200000]\n");
        printf("--window-step         INT      Step of sliding windows;                [10000]\n");
        printf("--min-marker          INT      Minimum number of markers of a window;      [2]\n");
        printf("--min-quality         INT      Minimum quality of markers for calculation [25]\n");
      //printf("--interval-max-mean   DOUBLE   Maximum mean of AFs;                     [1.00]\n");        
      //printf("--interval-min-mean   DOUBLE   Minimum mean of AFs;                     [0.99]\n");
      //printf("--interval-max-cvar   DOUBLE   Maximum CV of AFs;                       [0.01]\n");
        printf("-non-EMS                       Keep only EMS mutations;                  [off]\n");
        
        // options to plot proportion confidence interval at markers -- new on 2015-06-15 16:04
        printf("#Plot                          PCI in a defined mapping interval   \n");
        printf("-pci                           Turn on plot PCI;                         [off]\n");
        printf("--pci-cfd             DOUBLE   confidence level to plot PCI;            [0.95]\n");
        printf("--pci-chr             INT      ID of chr to plot PCI;                  [unset]\n");
        printf("--pci-start           INT      Position to start ploting PCI;          [unset]\n");
        printf("--pci-end             INT      Position to end ploting PCI;            [unset]\n");
        
        printf("-verbose                       Be talkative;                             [off]\n");
        printf("                                                                              \n");
        printf("*bg-background;fg-foreground; AF-Allele Frequency.                            \n");
        printf("*Filter 3 is more effective than filter 2.                                    \n");  
        printf("*PCI: proportion confidence interval.                                         \n");
        printf("*All SNPs are differences between (a) target genome(s) and a reference genome.\n");             
        printf("                                                                              \n");     
        exit(1);
    }
    
    /* set default initials */
    marker_score      =      25;
    marker_freq       =     0.0;
    reg_freq_min      =     0.2;   // for check_ref_err(...)
    filter_min_coverage =     1;   // == marker_read
    filter_max_coverage =   INF;
    marker_hit        =(double)INF;// markers with avg_hit > 1 are usually noises - but, caution!
    bg_score          =       1;
    bg_freq           =     0.0;
    bg_read           =       1;
    fg_N_cov          =     INF;   // <=> marker_N_cov;     INF means no filtering with this para
    fg_INDEL_cov      =     INF;   // <=> marker_INDEL_cov; INF means no filtering with this para
    bg_ref_base_file  =      "";
    bg_ref_cov        =      10;
    bg_ref_cov_max    =     500;
    bg_ref_freq       =     0.8;
    bg_ref_score      =      25;
    plot_bc           =       0;
    plot_boost        =       0;
    plot_marker       =       1;   // function for OC only
    plot_record       =       0;
    plot_scale        =       0;   // CMD["-plot-scale"]         = "";    
    only_EMS          =       0;
    window_size       =  200000;   // CMD["--window-size"] = "";
    window_step       =   10000;   // CMD["--window-step"] = "";
    interval_min_mean =    0.99;
    interval_max_mean =    1.00;
    interval_max_cvar =    0.01;
    quality_min       =    25.0;
    filter_min_marker =       2;   // CMD["--min-marker"]  = "";
    clusterK          =       5;   // number of clusters for visualization
    pci               = false;     // default: do not plot proportion confidence interval
    pci_chr           = "";        // chromosome id  to plot PCI
    pci_start         = 1;         // start position to plot PCI
    pci_end           = 1000;      // end   position to plot PCI
    pci_cfd           = 0.95;      // confidence level  for  PCI
    
    double startT= clock();
    bool   plot_monitor = true;
    /* 0.read practical values from cmd line, and write log accordingly */
    cmd_backcross(argc, argv);
    reg_freq_min      =   marker_freq;   // reset for check_ref_err(...)
    std::string fBClog("");
    fBClog = out_folder +  "BackCross.log";
    FILE* fplog = fopen(fBClog.c_str(), "w");
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
        fprintf(fplog, "Cannot open file to write log. Exited (in init_backcross(...)).\n");
        printf("Cannot open file to write log. Exited (in init_backcross(...)).\n"); exit(1);
    }
    fprintf(fplog, "\n\nOutputs achieved:\n");
    
    std::string fmarker2plot = "";                             // a  file  of   markers   to    plot
    /* 1.clean foreground markers according to background markers                                 */
    if(CMD.find("--bg") != CMD.end())
    {
        /* new 2014-03-18:  record  bg  cov  and  quality  for  future  clustering   of   markers */
        filter_fg_with_bg((char*)fbackground_file.c_str()); 
        /* frun_file will be "/path/to/SHOREmap_marker.bg_corrected"                              */
        fprintf(fplog, "\n\nFG markers in %s (filtered according to BG markers in %s):\n\t\t%s\n", 
                fmarker.c_str(), fbackground_file.c_str(), frun_file.c_str());
        if(verbose)
        printf("FG markers in %s (filtered according to BG markers in %s):\n\t\t%s\n", 
                fmarker.c_str(), fbackground_file.c_str(), frun_file.c_str());                
    }
    else
    {
        frun_file  = "";
        frun_file += fmarker;
        /* frun_file will be "/path/to/quality_variant.txt"                                       */
        fprintf(fplog, "\n\nForeground markers (not bg-corrected):\n\t\t%s\n", frun_file.c_str());
        if(verbose)
        printf("Foreground markers (not bg-corrected):\n\t\t%s\n", frun_file.c_str());        
    }
    if(plot_bc && bg_ref_filter==0) fmarker2plot = frun_file;
    
    /* 2.filter foreground markers: result in three filtered marker files                         */
    /*   *_q25, *_q25_f20, *_q25_f20_EMS will be generated: quality/ccd/ems filtered.             */
    std::string snp_1st_fout1 = "";
    std::string snp_1st_fout2 = "";
    std::string snp_1st_fout3 = "";
    int ft2 = filter_fg((char*)frun_file.c_str(), &snp_1st_fout1, &snp_1st_fout2, &snp_1st_fout3);
    if(ft2 == 1)
    {
        fprintf(fplog, "\nFG markers in %s ", frun_file.c_str());
        fprintf(fplog, "(filtered according to self-quality, concordance and EMS): \n");
        fprintf(fplog, "\t\t%s\n", snp_1st_fout1.c_str());
        fprintf(fplog, "\t\t%s\n", snp_1st_fout2.c_str());
        if(only_EMS==1)
        {
            fprintf(fplog, "\t\t%s\n", snp_1st_fout3.c_str());
        }
        
        if(plot_bc && bg_ref_filter==0 && only_EMS==1)
        {
            fmarker2plot = snp_1st_fout3;
            fprintf(fplog, "As EMS has been specified, visualized markers are in: %s\n. ", 
                    snp_1st_fout3.c_str());
        }
        else if(plot_bc && bg_ref_filter==0)
        {
            fmarker2plot = snp_1st_fout2;
            fprintf(fplog, "By default, if EMS is not specified, visualized markers are in: %s\n. ", 
                    snp_1st_fout2.c_str());
        }
    }
    else
    {
        printf("ERROR in 1st filtering (in init_backcross(...)). Exited.\n");
        exit(1);
    }
    /* 3.filter foreground markers produced in last step if cmd is set                            */
    if(bg_ref_filter == 1)
    {
        std::string fg_snp_file = "";
        if(ft2 == 1)
        {
            if(only_EMS == 1)
            {
                fg_snp_file  = snp_1st_fout3;
                fprintf(fplog, "As EMS has been specified, markers for ref-filter are in: %s\n. ", 
                    snp_1st_fout3.c_str());
            }
            else
            {
                fg_snp_file  = snp_1st_fout2;
                fprintf(fplog, "EMS not specified, markers for ref-filter are in: %s (default). \n", 
                    snp_1st_fout2.c_str());
            }
        }
        else
        {
            /* never happens: TO REMOVE */
            fg_snp_file      = frun_file;
        }
        if(verbose)
        {
            printf("Markers for ref-filter are in: %s.\n",  fg_snp_file.c_str());
        }
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
        fg_snp_map_filtered.clear();    // TODO - remove redundant = modify func to let it know when 
        fprintf(fplog,"\nFG markers in %s (filtered according to BG-ref-base quality, concordance):",
                fg_snp_file.c_str());
        fprintf(fplog, "\n\t\t%s\n", fg_snp_file_filtered.c_str());
        if(plot_bc)
        {
            fmarker2plot = fg_snp_file_filtered;
            fprintf(fplog, "\n\t\tMarkers in this file will be visualized. \n");
        }
    }
    /* 4. plot allele frequency if required                                                       */
    if(plot_bc)
    {
        /* 1. read markers                                                                        */
        unsigned long num_m_given = 0;
        unsigned long num_m_inuse = 0;
        if(!read_marker((char*)fmarker2plot.c_str(), &num_m_given, &num_m_inuse))
        {
            fprintf(fplog, "ERROR: no marker info recorded from file %s.\n", fmarker2plot.c_str());
            printf("ERROR: no marker info recorded. Exited. \n"); exit(1);
        }
        /* 2. read consensus base and allele1 allele2 error counts                                */
        if (!read_allele_counts2((char*)fconsensus.c_str()))
        {
            fprintf(fplog, "ERROR: no consensus info recorded. Program Exited. \n");
            printf("ERROR: no consensus info recorded. Exited. \n"); exit(1);
        }
        /* 3. visualize frequency of markers                                                      */
        printf("Start plotting allele frequency: \n\n");
        if(!visualizeFreq_bc())
        {
            fprintf(fplog,"\n\nVisualization failed on markers of %s. \n", fmarker2plot.c_str());
            printf("   ERROR: visualization failed. Exited. \n"); exit(1);
        }
        else
        {
            fprintf(fplog,"\n\nVisualization done on markers of %s. \n", fmarker2plot.c_str());
        }
    }
    
    /* 5.close log file                                                                           */
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    fprintf(fplog, "\n\nBackcross function successfully finished on %s\n", asctime(tinfo));
    double finishT= clock();
    fprintf(fplog, "\nTime consumed %.4f seconds.\n", (finishT-startT)/1000000);
    fclose(fplog);
    return 1;
}

/* visualize allele frequency of chromosomes one by one                                           */
bool visualizeFreq_bc()
{
    map<std::string, unsigned long> chrsize = CHR2SIZE;
    /* begin: visualize the chromosomes one by one                                                */
    /* set          up            pdf         file                                                */
    char chrNameBuf[512];
    if(plot_boost)
    sprintf(chrNameBuf, "%sBC_AF_visualization_%ld_boost_%s.pdf\0", out_folder.c_str(), runid, reg_chromosome.c_str());
    else
    sprintf(chrNameBuf, "%sBC_AF_visualization_%ld_%s.pdf\0",       out_folder.c_str(), runid, reg_chromosome.c_str());
    metafl("pdf");                         // level 0     - defines    the      metafile      format
    setfil(chrNameBuf);                    // level 0     - sets         alternate          filename
    filmod("VERSION");                     // level 0/1/2/3 do not shorten existing filename
    scrmod("revers");                      // level 0     - swaps  back   and   foreground   colours
    setpag("USEL");                        // level 0     - selects   a   predefined   page   format
                                           //               landscape:       11180*8640       points
    /* initialization                                                                             */
    disini();                              // level 0     - initialize DISLIN with default PARAs and
                                           //               a plotfile. The level  is  set   to   1. 
                                           //               DISINI must be called before  any  other 
                                           //               DISLIN   routine     except    for   ...
    unsigned long plot_ith_chr = 0;
    map<std::string, unsigned long>::iterator chr_itr;
    for(chr_itr = chrsize.begin(); chr_itr != chrsize.end(); chr_itr ++)
    {
        printf("   Visualizing allele frequency of chromosome %s\n", (*chr_itr).first.c_str());
        map<unsigned long, TRIPLE> internalData = CHR2POS2_ale1_ale2_err_COUNT[(*chr_itr).first];
        if(internalData.size() == 0) 
        {
            printf("No markers to visualize for chromosome %s.\n\n", (*chr_itr).first.c_str());
            continue;
        }
        /* begin: plot chromosome                                                                 */
        map<unsigned long, TRIPLE> filtered; /* not used                                          */
        filtered.clear();
        plot_ith_chr ++;
        if(plot_ith_chr > 1) newpag();
        plot_chr_winboost((*chr_itr).first, plot_ith_chr, internalData, filtered);
        printf("   Visualization of chromosome %s done. \n\n", (*chr_itr).first.c_str());
        /* end:   plot chromosome                                                                 */    
    }
    /* end:  visualizing the chromosomes one by one                                               */
    /* close pdf file */ 
    disfin();                             // level 1/2/3 - terminates DISLIN and prints a message on
                                          //               the screen. The level is set  back  to  0
                                          //               not    required      for      qplsca(...)    
    return true;
}

/* this function filters foreground markers according to background markers                       */
/* TODO: make this function capable of reading many background files. - 2013-03-25                */
int  filter_fg_with_bg(char* fileBG)
{    
    /* 1. read and record background marker info for future reference                             */
    FILE* FileBG = fopen(fileBG, "r");
    if(!FileBG)
    {
        printf("Background file \'%s\' does NOT exist. Exited.\n", fileBG);
        exit(1);
    }
    
    bgMARKER.clear();                      // this info will be used in future to cluster fg-markers
    unsigned long itemBG = 0;
    while(!feof(FileBG))
    {
        char projname[32];
        char chrid[32];
        char position[32];                                                     // caution, not long!
        char refb[32];
        char mutb[32];
        char qscore[32];
        char pcoverage[32];
        char pconcordance[32];
        char pavghits[32];          
        fscanf(FileBG, "%s\n", projname);
        fscanf(FileBG, "%s\n", chrid);
        fscanf(FileBG, "%s\n", position);
        fscanf(FileBG, "%s\n", refb);
        fscanf(FileBG, "%s\n", mutb);
        fscanf(FileBG, "%s\n", qscore);
        fscanf(FileBG, "%s\n", pcoverage);
        fscanf(FileBG, "%s\n", pconcordance);
        fscanf(FileBG, "%s\n", pavghits);
        
        if(mutb[0] == '-') continue;
        
        itemBG ++;
        //if( (atol(qscore)>=bg_score) && (atof(pcoverage)>=bg_read) && (atof(pconcordance)>=bg_freq))
        //{
        std::string stemp = (string)chrid+".#."+(string)position;
        std::string itemp = (string)qscore+"#"+(string)pcoverage+"#"+(string)pconcordance+"#"+(string)mutb;
        bgMARKER.insert(std::pair<string, string>(stemp, itemp));
        //}
    }
    fclose(FileBG);
    if(verbose)
    printf("%ld bg mutations recorded to filter fg markers.\n", itemBG);
    
    /* 2. select the foreground markers accroding to given background markers */
    char out_tem[1024];
    sprintf(out_tem, "%sSHOREmap_marker.bg_corrected\0", (char*)out_folder.c_str());
    FILE* fpcorrected = fopen(out_tem, "w");
    if(fpcorrected == NULL)
    {
        printf("cannot open file to write info: \t%s\n", out_tem);
        exit(1);
    }
    
    /* 2.1 check foreground markers */
    FILE* FileFG =  fopen(fmarker.c_str(), "r");
    if(!FileFG)
    {
        printf("Foreground file \'%s\' does NOT exist. Exited.\n", fmarker.c_str());
        exit(1);
    }
    
    if(verbose)
    printf("Traversing and filtering fg mutations...\n");
    bool firstw = true;
    std::string copied_line;
    while(!feof(FileFG))
    {
        char projname[32];
        char chrid[32];
        char position[32]; // caution, not long.
        char refb[32];
        char mutb[32];
        char qscore[32];
        char pcoverage[32];
        char pconcordance[32];
        char pavghits[32];
              
        fscanf(FileFG, "%s\n", projname);
        fscanf(FileFG, "%s\n", chrid);
        fscanf(FileFG, "%s\n", position);
        fscanf(FileFG, "%s\n", refb);
        fscanf(FileFG, "%s\n", mutb);
        fscanf(FileFG, "%s\n", qscore);
        fscanf(FileFG, "%s\n", pcoverage);
        fscanf(FileFG, "%s\n", pconcordance);
        fscanf(FileFG, "%s\n", pavghits);
        
        /* ignore chrosomome not in list */
        if(CHR2SIZE.find((std::string)chrid) == CHR2SIZE.end()) continue;
        /* ignore deletion variation '-' */
        if(mutb[0] == '-') continue;
        
        string  bgkey("");
        bgkey   = (string)chrid+".#."+(string)position;
        map<string, string>::iterator bgmkr_itr;
        bgmkr_itr = bgMARKER.find(bgkey);
        
        bool myfilter = false;
        if(bgmkr_itr != bgMARKER.end())
        {
            string bgmkrtmp = (*bgmkr_itr).second;
            std::vector<string> bgmkr_info = split_string(bgmkrtmp, '#');
            //if( (atol(qscore)>=bg_score) && (atof(pcoverage)>=bg_read) && (atof(pconcordance)>=bg_freq))
            unsigned long bgscore_tmp = atol(bgmkr_info[0].c_str());
            double         bgread_tmp = atof(bgmkr_info[1].c_str());
            double         bgccod_tmp = atof(bgmkr_info[2].c_str());
            if(bgmkr_info[3].find(string(mutb))==0) // when two mut bases are the same, remove mkr
            {
                if((bgscore_tmp>=bg_score) && (bgread_tmp>=bg_read) && (bgccod_tmp>=bg_freq))
                {
                    /* need to remove the corresponding fg marker */
                    myfilter = true;
                    /* remove this bg info     */
                    bgMARKER.erase(bgmkr_itr);
                } 
            }
            else
            {
                ;
                //printf("same position, diff mutb: %s\t%s\t%s\t%c\t%c\t%s\t%s\t%s\t%s\n", 
                //  projname, chrid, position, refb[0], mutb[0], qscore, pcoverage, pconcordance, pavghits);
            }
        }
        
        if(bgmkr_itr == bgMARKER.end() || !myfilter)
        {
            if(atof(pavghits)<=marker_hit) // caution: filtering with avg_hit of a candidate snp
            {
                if(firstw == true)
                {
                    firstw = false;
                }
                else
                {
                    fprintf(fpcorrected, "\n");
                }
                fprintf(fpcorrected, "%s\t%s\t%s\t%c\t%c\t%s\t%s\t%s\t%s", 
                  projname, chrid, position, refb[0], mutb[0], qscore, pcoverage, pconcordance, pavghits);
                
                /* set bg marker with information on mut base as 0+#+0+#0 */
                if(bgmkr_itr == bgMARKER.end())
                {
                    bgMARKER.insert(std::pair<string, string>(bgkey, "0#0#0"));
                }
            }
        }
    }
    fclose(FileFG);
    fclose(fpcorrected);
    printf("Traversing and filtering fg mutations done.\n");
    
    /* set path and filename of markers to be used for filtering */
    frun_file  = "";
    frun_file += out_tem;                                 // out_folder/SHOREmap_marker.bg_corrected
    
    return 1;
}

/* initialize parameters according to the cmd line */
int  cmd_backcross(int argc, char* argv[])
{
    int ic = 2;
    while (ic < argc)
    {
        if(!strcmp(argv[ic],"-verbose"))                  // option to variable: verbose
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            verbose         = 1;
            printf("Be talkative during process. \n");
        }
        ic ++;
    }
    
    ic = 2;                                               // option ic=0: SHOREmap; ic=1: backcross
    while (ic < argc) 
    {
        strcatCMD += " ";
        strcatCMD += (string)argv[ic];
        if(!strcmp(argv[ic],"--chrsizes"))                // option 1 to variable: fchrsizes
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fchrsizes         = argv[ic];                 // set value for option
            FILE* fp_chrsizes = fopen((char*)fchrsizes.c_str(), "r");
            if(!fp_chrsizes)  { print_error_exit(argv[ic-1], false); }
            fclose(fp_chrsizes);                          // check done.
            if (verbose) printf("Chrs sizes read from file:\t%s\n", fchrsizes.c_str());
            if(!read_chromosomes((char*)fchrsizes.c_str()))// read contents
            {
                printf("ERROR: invalid content in file %s.\n", fchrsizes.c_str());
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--folder"))             // option 2 to variable: out_folder
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            out_folder          = argv[ic];
            if(out_folder[out_folder.length()-1] != '/') out_folder += "/";
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
                    printf("ERROR: cannot create output path. Exited. %s\n", out_folder.c_str());
                    exit(1);
                }
            }
            closedir(dir_out_folder);
        }
        else if(!strcmp(argv[ic],"--marker"))              // option 3 to variable: fmarker
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker         = argv[ic];
            FILE* fp_marker = fopen((char*)fmarker.c_str(), "r"); 
            if(fp_marker == NULL)
            {
                printf("marker file \'%s\' does NOT exist. Exited.\n", fmarker.c_str());
                exit(1);
            }
            if (verbose) printf("File of markers provided:\t%s\n", fmarker.c_str());
           // read_marker();
            fclose(fp_marker);
        }
        else if(!strcmp(argv[ic],"--marker-score"))       // option 4 to variable: marker_score
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            marker_score           = atol(argv[ic]);
            if(marker_score < 1) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("marker_score:\t\t%ld\n", marker_score);
        }
        else if(!strcmp(argv[ic],"--marker-freq"))        // option 5 to variable: marker_freq
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            marker_freq     = atof(argv[ic]);
            if(marker_freq<0.0 || marker_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("marker_freq:\t\t%.2f.\n", marker_freq);
        }
       else if(!strcmp(argv[ic],"--min-coverage"))        // option 6 tovariable:filter_min_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_min_coverage set as:\t\t%ld\n", filter_min_coverage);
        }
        else if(!strcmp(argv[ic],"--max-coverage"))       // option 6.1tovariable:filter_max_coverage
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_max_coverage   = atol(argv[ic]);
            if (verbose) printf("filter_max_coverage set as:\t\t%ld\n", filter_max_coverage);
        }
        else if(!strcmp(argv[ic],"--marker-hit"))         // option 6.5 to variable: marker_hit
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            marker_hit           = atof(argv[ic]);
            if(marker_hit<0) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("marker_hit:\t\t%.4f\n", marker_hit);
        }
        else if(!strcmp(argv[ic],"--bg"))                 // option 7 to variable: fbackground_file
        {
            // TODO: many background files might be provided. - 2013-04-16 21:50
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fbackground_file = argv[ic];
            FILE* fp_bg      = fopen((char*)fbackground_file.c_str(), "r"); 
            if(fp_bg        == NULL)
            {
                printf("background file \'%s\' does NOT exist. Exited.\n", fbackground_file.c_str());
                exit(1);
            }
            if (verbose) printf("File of background markers:\t%s\n", fbackground_file.c_str());
            // read_bg_marker();
            fclose(fp_bg);
        }
        else if(!strcmp(argv[ic],"--bg-score"))           // option 8 to variable: bg_score
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_score    = atol(argv[ic]);
            if(bg_score<1) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("bg_score   :\t\t%ld\n", bg_score);
        }
        else if(!strcmp(argv[ic],"--bg-freq"))            // option 9 to variable: bg_freq
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_freq     = atof(argv[ic]);
            if(bg_freq<0.0 || bg_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("bg_freq    :\t\t%.2f.\n", bg_freq);
        }
        else if(!strcmp(argv[ic],"--bg-cov"))             // option 10 to variable: bg_read
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_read           = atol(argv[ic]);
            if(bg_read<1) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("bg_read    :\t\t%ld\n", bg_read);
        }
        else if(!strcmp(argv[ic],"-no-summary"))          // option 11 to variable: summmary
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            summary       = 0;
            if (verbose) printf("summary    :\t\t%ld\n", summary);
        }
       else if(!strcmp(argv[ic],"-non-EMS"))              // option 12 to variable: only_EMS
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            only_EMS = 1;
            if (verbose) printf("only_EMS   :\t\t%ld\n", only_EMS);
        }
        else if(!strcmp(argv[ic],"-no-filter"))           // option 13 to variable: filter_plot
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            filter_plot = 0;
            if (verbose) printf("filter_plot:\t\t%ld\n", filter_plot);
        }
        else if(!strcmp(argv[ic],"-other-mutant"))        // option 14 to variable: other_mutant
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            other_mutant = 0;
            if (verbose) printf("other_mutant set as:\t\t%ld\n", other_mutant);
        }
        else if(!strcmp(argv[ic],"-bg-ref-filter"))       // option 14.1 to variable: bg_ref_filter
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            bg_ref_filter = 1;
            if (verbose) printf("bg_ref_filter set as:\t\t%ld\n", bg_ref_filter);
        }
        else if(!strcmp(argv[ic],"-plot-bc"))             // option 14.2 to variable: plot_bc
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_bc = 1;
            if (verbose) printf("plot_bc set as:\t\t%ld\n", plot_bc);
        }
        else if(!strcmp(argv[ic],"-plot-scale"))          // option 14.25 to variable: plot_scale.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_scale     = 1; // caution
            if (verbose) printf("plot_scale set as:\t\t\ttrue\n", plot_scale);
        }
        else if(!strcmp(argv[ic],"--window-size"))        // option 14.3 to variable: window_size
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_size           = atol(argv[ic]);
            if(window_size<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("window_size set as:\t\t\t%ld\n", window_size);
        }
        else if(!strcmp(argv[ic],"--window-step"))        // option 14.4 to variable: window_step
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            window_step           = atol(argv[ic]);
            if(window_step<=1) 
                { printf("ERROR: arg of %s must be larger than 1. Exited.\n",
                       argv[ic-1]); exit(1);}
            if (verbose) printf("window_step set as:\t\t\t%ld\n", window_step);
        }
        else if(!strcmp(argv[ic],"--min-marker"))      // option 14.5 to variable: filter_min_marker
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            filter_min_marker   = atol(argv[ic]);
            if(filter_min_marker<1) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("filter_min_marker set as:\t\t%ld\n", filter_min_marker);
        }
        else if(!strcmp(argv[ic],"--consen"))             // option 15.0 to variable: fconsensus
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fconsensus         = argv[ic];
            FILE* fp_consensus = fopen((char*)fconsensus.c_str(), "r"); 
            if(fp_consensus   == NULL)
            {
                printf("consensus file \'%s\' does NOT exist. Exited.\n", fconsensus.c_str());
                exit(1);
            }
            fclose(fp_consensus);
        }/*
        else if(!strcmp(argv[ic],"--bg-ref-consen"))      // option 15.0 to variable: fconsenref
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fconsenref         = argv[ic];
            FILE* fp_consenref = fopen((char*)fconsenref.c_str(), "r"); 
            if(fp_consenref   == NULL)
            {
                printf("ref consensus file \'%s\' does NOT exist. Exited.\n", fconsenref.c_str());
                exit(1);
            }
            fclose(fp_consenref);
        }*/
        else if(!strcmp(argv[ic],"--fg-N-cov"))          // option 16 to variable: fg_N_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            fg_N_cov   = atol(argv[ic]);
            if (verbose) printf("fg_N_cov set as:\t\t%ld\n", fg_N_cov);
        }
        else if(!strcmp(argv[ic],"--fg-INDEL-cov"))      // option 17 to variable: fg_INDEL_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            fg_INDEL_cov   = atol(argv[ic]);
            if (verbose) printf("fg_INDEL_cov set as:\t\t%ld\n", fg_INDEL_cov);
        }
        else if(!strcmp(argv[ic],"--bg-ref-base-file"))  // option 18 to variable: bg_ref_base_file
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
        else if(!strcmp(argv[ic],"--bg-ref-cov"))        // option 19 to variable: bg_ref_cov
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov set as:\t\t%ld\n", bg_ref_cov);
        }
        else if(!strcmp(argv[ic],"--bg-ref-cov-max"))    // option 35.5 to variable: bg_ref_cov_max
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_cov_max   = atol(argv[ic]);
            if (verbose) printf("bg_ref_cov_max set as:\t\t%ld\n", bg_ref_cov_max);
        }
        else if(!strcmp(argv[ic],"--bg-ref-freq"))       // option 20 to variable: bg_ref_freq
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_freq     = atof(argv[ic]);
            if(bg_ref_freq<0.0 || bg_ref_freq>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1); }
            if (verbose) printf("bg_ref_freq:\t\t%.2f.\n", bg_ref_freq);
        }
        else if(!strcmp(argv[ic],"--bg-ref-score"))      // option 21 to variable: bg_ref_score
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            bg_ref_score           = atol(argv[ic]);
            if(bg_ref_score<1) 
                { printf("ERROR: arg of %s must be larger than 0. Exited.\n",
                      argv[ic-1]); exit(1);}
            if (verbose) printf("bg_ref_score:\t\t%ld\n", bg_ref_score);
        }
        else if(!strcmp(argv[ic],"-plot-boost"))         // option 22 to variable: plot_boost
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_boost     = 1;                          // caution
            if (verbose) printf("plot_boost set as:\t\t%d\n", plot_boost);
        }
        else if(!strcmp(argv[ic],"-plot-win"))           // option 22.5 to variable: plot_window
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_window     = 1;                         // caution
            if (verbose) printf("plot_window set as:\t\t%d\n", plot_window);
        }
        else if(!strcmp(argv[ic],"-rab"))                // option 38 to variable: plot_record.......
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            plot_record     = 1;
            if (verbose) printf("plot_record set as:\t\t\t%ld\n", plot_record);
        }
        else if(!strcmp(argv[ic],"-runid"))              // option 23 to variable: runid
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            runid = atoi(argv[ic]);
            if (verbose) printf("runid has set as:\t\t\t%ld.\n", runid);
        }
        else if(!strcmp(argv[ic],"--interval-max-mean")) // option 23.5 tovariable:interval_max_mean
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_max_mean       = atof(argv[ic]);
            if (verbose) 
            {
               printf("interval_max_mean set as:\t\t%.2f by option %s.\n", 
                      interval_max_mean, argv[ic-1]);
            }
        }
        else if(!strcmp(argv[ic],"--interval-min-mean")) // option 24 to variable: interval_min_mean
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_min_mean       = atof(argv[ic]);
            if (verbose) 
            {
               printf("interval_min_mean set as:\t\t%.2f by option %s.\n", 
                      interval_min_mean, argv[ic-1]);
            }
        }
        else if(!strcmp(argv[ic],"--interval-max-cvar")) // option 25 to variable: interval_max_cvar
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            interval_max_cvar    = atof(argv[ic]);
            if (verbose)
            {
                printf("interval_max_cvar set as:\t\t%.2f by option %s.\n", 
                        interval_max_cvar, argv[ic-1]);
            }
        }
        else if(!strcmp(argv[ic],"--min-quality"))// option 26 to variable: quality_min
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            quality_min    = atof(argv[ic]);
            if (verbose)
            {
                printf("quality_min set as:\t\t%.2f by option %s.\n", 
                        quality_min, argv[ic-1]);
            }
        } 
        else if(!strcmp(argv[ic],"--cluster"))    // option 27 to variable: clusterK
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            clusterK   = atol(argv[ic]);
            if(clusterK<1 || clusterK>40) 
                { printf("ERROR: arg of %s must be in [1, 40]. Exited.\n", 
                       argv[ic-1]); exit(1);}
            if (verbose) printf("clusterK set as:\t\t%ld\n", clusterK);
        }
        else if(!strcmp(argv[ic],"-pci"))        // option 28 to variable: pci
        { 
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            pci = true;
            if (verbose) printf("pci set as:\t\t%ld\n", pci);
        }
        else if(!strcmp(argv[ic], "--pci-chr"))  // option 29 to variable: pci_chr....
        {
            ic       = precheck_opt(argc, argv, ic, argv[ic], true, false);
            pci_chr += argv[ic];
            if (verbose) printf("pci_chr set as:\t\t%s.\n", pci_chr.c_str());
            if(CHR2SIZE.find(pci_chr)==CHR2SIZE.end()) 
            {
                printf("ERROR: chromosome %s to plot PCI is not found in the list. Exited.\n", 
                    (char*)pci_chr.c_str()); 
                exit(1);
            }
        }
        else if(!strcmp(argv[ic],"--pci-start")) // option 30 to variable: pci_start.........
        {
            if(pci_chr == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pci_start      = atol(argv[ic]);
            string schr(pci_chr);
            if(pci_start<0 || pci_start>CHR2SIZE[schr])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                     argv[ic-1]); exit(1);}
            if (verbose) printf("pci_start set as:\t\t\t%ld\n", pci_start);
        }
        else if(!strcmp(argv[ic],"--pci-end")) // option 31 to variable: pci_end.........
        {
            if(pci_chr == "") {
                printf("Chromosome ID has not been set.\n"); exit(1);}
            
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pci_end      = atol(argv[ic]);
            string schr(pci_chr);
            if(pci_end<0 || pci_end>CHR2SIZE[schr])
                { printf("ERROR: arg of %s exceeds bounds of chromosome. Exited.\n", 
                     argv[ic-1]); exit(1);}
            if (verbose) printf("pci_end set as:\t\t\t%ld\n", pci_start);
        }
        else if(!strcmp(argv[ic],"--pci-cfd")) // option 32 to variable: pci_cfd
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            pci_cfd    = atof(argv[ic]);
            if (verbose)
            {
                printf("pci_cfd set as:\t\t%.2f.\n", pci_cfd);
            }
        }
        else if(!strcmp(argv[ic],"-verbose"))     //
        {
            ;
        }
        else
        {
            printf("Warning: parameter NOT recognized: %s.\n", argv[ic]);
        }
        ic ++;
    }
    /* Check: mandatory parameters; or optional parameters should be provided simultaneously      */
    bool mandatoryTF  = true;
    bool mandatoryTF2 = true;
    if(CMD.find("--chrsizes")    == CMD.end())      {mandatoryTF = false;}
    else if(CMD.find("--folder") == CMD.end())      {mandatoryTF = false;}
    else if(CMD.find("--marker") == CMD.end())      {mandatoryTF = false;}
    else;
    
    if(bg_ref_filter)
    if(CMD.find("--bg-ref-base-file") == CMD.end()) {mandatoryTF2 = false;}
    if(bg_ref_filter)
    if(CMD.find("--consen")    == CMD.end())        {mandatoryTF2 = false;}  
    if(CMD.find("-plot-bc")    != CMD.end()) 
    if(CMD.find("--consen")    == CMD.end())        {mandatoryTF2 = false;}
    
    if(!mandatoryTF)
    {
        printf("Parameters required for --chrsizes/folder/marker. Exited. (in init_backcross(...))\n");
        exit(1);
    }
    if(!mandatoryTF2)
    {
        printf("Consen info file of markers required. Exited. (in init_backcross.cpp) \n");
        exit(1);
    }
    
    if(filter_min_coverage > filter_max_coverage)
    {
        printf("ERROR: min and max coverages for filtering markers do not agree with each other.\n");
        exit(1);
    }
    
    cluster_avg_coverage = (filter_min_coverage + filter_max_coverage)/2;
    cluster_dim_coverage = (filter_max_coverage - filter_min_coverage)/2;
    if(cluster_dim_coverage == 0) cluster_dim_coverage = 1; // caution
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
    
    if(other_mutant==1 && only_EMS==1) only_EMS = 0;
    return 1;
}
