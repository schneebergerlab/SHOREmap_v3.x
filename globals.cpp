/* globals.h: declare the global variables that are used through the whole program. ****************************************************/
#include      <map>                                                                                                                    //
#include   <string>                                                                                                                    //
#include <stddef.h>                                                                                                                    //
                                                                                                                                       //
#define INF      9999999999
#define MARGIN 0.0000000001                                                                                                                                       
using namespace std;                                                                                                                   //
struct TRIPLE
{
    unsigned long Ci[3];
};
struct WINSAVG
{
    unsigned long wbegin;
    unsigned long wend;
    double wfreq;
    double wr;                       // r-value not used
    double wboost;
};
struct LIKELIHOOD                    // not-used
{
    double P1est;
    double errEst;
    double min;
};
struct INTERVAL
{
    unsigned long start;
    unsigned long size;
    double value;
    double lev;                      // level
};
struct REGION                        // for checking storage_shoremapmle.
{
    unsigned long start;
    unsigned long size;
};

struct MARK6                         // in check_ref_err
{
    std::string ref;
    std::string mut;
    std::string qua;
    std::string cov;
    std::string ccd;
    std::string aht;
    std::string par;                 // indicates it is determined in which parent
};
/**************************************************************** Command line parameters **********************************************/
unsigned long row_first;             /* the first line to extract from the consensus file               parameters                     */ 
unsigned long row_last;              /* the last  line to extract from the consensus file               for                            */
unsigned long col_first;             /* the first row  to extract from the consensus file               extracting                     */
unsigned long col_last;              /* the last  row  to extract from the consensus file               consensus info                 */
unsigned long extract_bg_ref=0;      /* extracting from quality_reference.txt of a sequenced background mutant                         */
/***************************************************************************************************************************************/
char*         mycross;               /* outcross or backcross                                                                          */
double        expect = 1.0;          /* expected coverage at a locus - not used                                                        */
double        confidence;            /* confidence level [default 0.99].					                       */	
double      interval_min_mean = 0.99;/* minimum mean                     of allele frequency of a set of markers in a mapping interval */
double      interval_max_mean = 1.00;/* maximum mean                     of allele frequency of a set of markers in a mapping interval */
double      interval_max_cvar = 0.01;/* maximum coefficient of variation of allele frequency of a set of markers in a mapping interval */
std::string   fchrsizes   = "";      /* name of a file   which contains chromosome id \t size; to be initialized with input parameters.*/
std::string   out_folder  = "";      /* name of a folder where output files can be found; to be initialized wiht input parameters.     */	
std::string   fmarker     = "";      /* name of a file   which contains markers							       */
std::string   fmarker_pa  = "";      /* name of a file   which contains markers of parent a (used for filtering - not necessary)       */
std::string   fmarker_pb  = "";      /* name of a file   which contains markers of parent b (used for filtering - not necessary)       */
std::string   fconsensus  = "";      /* name of a file   which contains alleles							       */	
std::string   fconsenref  = "";      /* name of a file   which contains alleles	for a reference used for filtering                     */	
std::string   fpeaks      = "";      /* name of a file   which contains peak-positions after calling visualize or backcross function   */
unsigned long chrsizes_max;          /* maximum length of chr sizes                                                                    */ 
double        quality_min = 25.0;    /* minimum base quality score (for backcross only)                                                */  
double        quality_max = 40.0;    /* maximum base quality score used in find_marker_pos_2parents_v3.cpp                             */
unsigned long window_size;           /* window size for visualization and peak finding etc [default 50,000].		               */	
unsigned long window_step;           /* window step for ...[default 10,000].							       */
                                     /*                                                                                                */
unsigned long peak_window_size;      /* window size for peaking finding [default 50,000].					       */	
unsigned long peak_window_step;      /* window step for ...             [default 10,000].					       */
                                     /*                                                                                                */
unsigned long filter_min_marker;     /* threshold value: if a window has less markers than this minimum, it will be filtered.	       */
unsigned long filter_min_coverage;   /* threshold value: if a locus has a coverage less than this minimum, it will be filtered.        */	
unsigned long filter_max_coverage;   /* threshold value: if a locus has a coverage more than this maximum, it will be filtered.	       */
unsigned long cluster_avg_coverage;  /* threshold value: this value will be used for clustering of markers before plotting             */
unsigned long cluster_dim_coverage;  /* threshold value: this value will be used for clustering of markers before plotting as diameter */
                                     /*                                                                                                */
unsigned long outlier_window_size;   /* window size to assess local allel frequency used for outlier removal [default 200,000].        */
double        outlier_pvalue;        /* p_value for outlier removal.								       */	
                                     /*                                                                                                */
double        misphenotyped;         /* degree of putatively mis-scored plants [default 0.00]. ???			               */	
                                     /*                                                                                                */
std::string   reg_chromosome="";     /* name/id of chromosome to annotate/zoom. [caution: SHOREmap.pl cannot recognize string-ids].    */
unsigned long reg_begin;             /* to annotate/zoom from this point.						               */	
unsigned long reg_end;               /* to annotatezoom   to this point.							       */	
double        reg_freq_min;          /* minimum frequency of zoom region.							       */
double        reg_freq_max;          /* maximum frequency of zoom region.						               */
                                     /*                                                                                                */ 
std::string   freferror   = "";      /* name of a file which contains reference errors.	= parent-b              		       */	
unsigned long background2 = 0;       /* the mutation is in the second parent if it equals 1 [default 0].        		       */	
unsigned long verbose;               /* to be talktive or not during process.							       */
                                     /*                                                                                                */
unsigned long runid       = 1;       /* running id of a specific task.						        	       */
unsigned long r_max;                 /* distance of a specific mutation from a peak.						       */
bool          plot_r      = 0;       /* to plot frequency calculation 'r' instead of 'boost'. - to removed         		       */	
double        boost_max   = 0;       /* maximum boost value.									       */
bool          plot_boost  = 1;       /* to plot 'boost'.							     		       */	
bool          plot_marker = 1;       /* to add single markers in the plot. (caution! function for OC only)	       		       */
unsigned long plot_bc;               /* switch on or off plotting of markers from backcross process                                    */
bool          plot_scale  = false;   /* plot the chromosomes scaled to the one with the maximum length                                 */
bool          plot_window = false;   /* plot window-averaged allele frequency of markers                                               */
unsigned long plot_record = 0;       /* record statsitcal info that has been plotted in pdf                                            */
                                     /*                                                                                                */							
unsigned long marker_score = 25;     /* minimum score cutoff for filtering targeted (F2) markers [default 25].                         */
//unsigned long marker_read  = 0;    /* minimum read support [default 0].                                                              */    
double        marker_freq  = 0.2;    /* minimum concordance: ratio of reads supporting a feature to total coverage [default 0.20].     */
double        marker_hit   = (double)INF;/* maximum average number of hit about a marker                                               */  
unsigned long pmarker_score    = 0;  /* minimum score cutoff of parent-SNPs for filtering F2-markers [default 0-no par-filtering].     */
unsigned long pmarker_min_cov  = 0;  /* minimum coverage cutoff of parent-SNPS for filtering F2-markers [default 0-no par-filtering].  */
unsigned long pmarker_max_cov  = INF;/* maximum coverage cutoff of parent-SNPs for filtering F2-markers [default INF-no par-filtering].*/
double        pmarker_min_freq = 0;  /* minimum frequency of parent-SNPs for filtering F2-markers [default 0.0-no par-filtering]       */
std::string   pmarker_ab_ratio;      /* ratio of SNPs of parents [default 1,1,1,1 - the cutoffs for parent-SNPs are the same]          */          
std::string   frun_file    = "";     /* to indicate which F2-marker file is in use.			         		       */
unsigned long indel_size   = 0;      /* size of indels to extract from a vcf file 						       */

/* filtering according to ref info   */
unsigned long bg_ref_filter    = 0;  // trun on this filter or not
unsigned long fg_N_cov         = INF;// if not filtering with this para
unsigned long fg_INDEL_cov     = INF;// if not filtering with this para
std::string   bg_ref_base_file =  "";
std::string   bg_ref_base_file_pb = "";
unsigned long bg_ref_cov          =  10;
unsigned long bg_ref_cov_max      = 500;
double        bg_ref_freq         = 0.8;
unsigned long bg_ref_score        =  25;
/* filtering according to background SNPs info                                                                                         */		
std::string   fbackground_file ="";   /* name of files which contain backgroud mutations (TODO).                                       */	
unsigned long bg_score         = 0;   /* minimum score cutoff for filtering background SNPs [default 0].			       */	
unsigned long bg_read          = 0;   /* minimum read support for background SNPs [default 0].	          			       */	
double        bg_freq          = 0.20;/* minimum concordance for background SNPs [default 20].		         		       */	
/*                                                                                                                                     */
/********* ploting related variables ***************************************************************************************************/
/*                                                                                                                                     */
unsigned long clusterK     = 5;      /* levels of scores for ranking markers in visualization                                          */
unsigned long summary      = 1;      /* turn on ploting all chromosomes in a single page as summary.			               */	
unsigned long only_EMS     = 0;      /* if 1, only plot EMS-induced mutations 					        	       */	
unsigned long other_mutant = 0;      /* other mutagen is not used for screening					         	       */
unsigned long filter_plot  = 1;      /* to plot after filtering						         		       */
std::string   R_input = "";          /* [caution: this is removed in this version using C].					       */	
std::string   R_out   = "";          /* folder of output when ...								       */
int           vis_annotation = 0;    /* visualization of annotation along chromosomes                                                  */
/*                                                                                                                                     */
/******** Additional global variables **************************************************************************************************/
/*                                                                                                                                     */
map<std::string, unsigned long> CHR2SIZE;/* where chromosome and its size is recorded: CHR2SIZE["chr"]     = (unsigned long)size.      */
map<std::string, bool> REFERROR;         /* this marker is an error:                   REFERROR["chr.#.pos"] = 1. - not used           */
/* TODO: combine ALLELE1, ALLELE2, QUALITY1 as one map<std::string, std::vector<std::string> >                                         */	
map<std::string, std::string> ALLELE1;   /* where the mutated allele is recorded:      ALLELE1["chr.#.pos"]  = 'A/C/G/T'.	       */	
map<std::string, std::string> ALLELE2;   /* where the ref allele is recorded:          ALLELE2["chr.#.pos"]  = 'A/C/G/T'.	       */
map<std::string, std::string> QUALITY1;  /* where quality of base call of the mutated allele is recorded: QUALITY1["chr.#.pos"] =      */	
map<std::string, std::string> CMD ;      /* where the command line parameters are recorded.  					       */
std::string  strcatCMD  = "SHOREmap ";   /* a string concatenation of the cmds: TODO                                                   */	
std::string  pwd;                        /* path of working folder.								       */
map<std::string, std::string> bgMARKER;  /* information about bg-markers which have been used to filtering/retaining foreground markers*/
map<std::string, std::string> bgREF;     /* information about bg reference bases which have been used to keep a fg-marker              */
map<std::string, std::string> bgREFB4A;  /* information about bg ref bases of parent B used to keep a fg-marker in parent A - outcross */
map<std::string, std::string> bgREFA4B;  /* information about bg ref bases of parent A used to keep a fg-marker in parent B - outcross */
bool keep_common = false;                /* whether keep common snps of parental lines in marker list - create                         */
/*                                                                                                                                     */
/***** Collect marker counts and sliding window values *********************************************************************************************************************/
map<std::string, map<std::string, unsigned long> > CHR2POS2ALLELE1_COUNT; /* where the count for allele 1 is recorded: CHR2POS2ALLELE1_COUNT["chr"]["pos"] = allele1_count.*/	
map<std::string, map<std::string, unsigned long> > CHR2POS2ALLELE2_COUNT; /* where the count for allele 2 is recorded: CHR2POS2ALLELE2_COUNT["chr"]["pos"] = alelle2_count.*/
map<std::string, map<std::string, unsigned long> > CHR2POS2ERROR_COUNT;   /* where the count for error    is recorded: CHR2POS2ERROR_COUNT["chr"]["pos"]   = error_count.  */
map<std::string, map<unsigned long, TRIPLE> >      CHR2POS2_ale1_ale2_err_COUNT; /* shq: postion records like allele 1, allele 2, error - 2013-03-10 19:31  		   */
/* file for markers annotation */
std::string snp_file    = "";
std::string ins_file    = "";
std::string del_file    = "";
std::string gff_file    = "";
std::string genome_file = "";
/* info about a single SNP     */
struct mySNP
{
    std::string   ecotype;      //1
    std::string   chromosome;   //
    unsigned long position;     //
    std::string   stype;        //
    std::string   gene_id;      //5
    unsigned long gene_pos;     // 
    unsigned long cds_pos;      //
    unsigned long codon_pos;    //
    unsigned long ns_change;    //
    unsigned long new_stop;     //10
    unsigned long lost_stop;    //
    unsigned long splicechange; //
    std::string   ref_base;     //
    std::string   new_base;     //
    std::string   ref_aa;       //15
    std::string   new_aa;       //
    //map<std::string, std::string> domain_change;
    unsigned long support;      //
    double        concordance;  //
    double        quality;      //
    unsigned long peak_distance;//20
    double        marker_ratio; //
    std::string   source;       //22
};
multimap<unsigned long, mySNP> SNPlist; // SNP list to annotat with gene info (gff file)

/* info about an insertion or deletion */
struct myINDEL
{
    std::string   ecotype;
    std::string   chromosome;
    unsigned long begin;
    unsigned long end;
    std::string   seq;
    std::string   stype;
    std::string   gene_id;
    unsigned long gene_pos;  // 
    unsigned long cds_pos;   // 
    unsigned long codon_pos; //
    unsigned long ns_change;
    unsigned long new_stop;
    unsigned long lost_stop;
    unsigned long splicechange;
    //map<std::string, std::string> domain_change;
    unsigned long support;
    double        concordance;
    double        quality;
    unsigned long peak_distance;
    double        marker_ratio;
    std::string   source;
};
multimap<unsigned long, myINDEL> INSlist; // INSERTIONS to annotate with gene info (gff file)
multimap<unsigned long, myINDEL> DELlist; // DELETIONS  to annotate with gene info (gff file)

/* info about gene annotation */
struct QUARTET                            // it is not quartet anymore as frame info is added -- 2016-MAR-03.
{
    unsigned long start;
    unsigned long end;
    std::string   orien;
    std::string   seq;
    int           frame;
};
struct RANGE
{
    unsigned long start;
    unsigned long end;
};
map<std::string, QUARTET> gene_ann;
map<std::string, map<unsigned long, QUARTET> > coding_ann;
multimap<std::string, multimap<std::string, RANGE> > seq_type;
int first_codon_frame = 0;                                // frame value for the first codon: can be '0','1','2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
int last_codon_frame  = 0;                                // frame value for the last codon: can be '0','1','2'.

/* info about GeneSNPlist */
struct GeneSNPlist
{
    std::string gene_id;                                  //1
    std::string isoform;                                  // 
    std::string chromosome;                               //
    unsigned long start;                                  //
    unsigned long end;                                    //5                
    unsigned long cds_length;                             //           
    unsigned long protein_length;                         //
    std::string   orientation;                            //
    multimap<unsigned long, mySNP> SNP_lists;             //
    map<unsigned long, unsigned long> coding_SNP;         //10
    map<unsigned long, unsigned long> NS_changes;         //
    map<unsigned long, unsigned long> AA_changes;         //
    map<unsigned long, unsigned long> CDSexon;            //
    map<std::string, std::string> protein;                //
    std::string ref_gene;                                 //15
    std::string ref_coding;                               //
    std::string eco_coding;                               //17
};

/* log with hints for solving possible errors */
std::string errLog = ".errlog"; // TODO 2013-04-19 09:56

/* z values of confidence interval @http://en.wikipedia.org/wiki/Standard_normal_table:Cumulative */
double z_table[32][11] =
{
    {0.0,0.50000,0.50399,0.50798,0.51197,0.51595,0.51994,0.52392,0.52790,0.53188,0.53586},
    {0.1,0.53980,0.54380,0.54776,0.55172,0.55567,0.55966,0.56360,0.56749,0.57142,0.57535},
    {0.2,0.57930,0.58317,0.58706,0.59095,0.59483,0.59871,0.60257,0.60642,0.61026,0.61409},
    {0.3,0.61791,0.62172,0.62552,0.62930,0.63307,0.63683,0.64058,0.64431,0.64803,0.65173},
    {0.4,0.65542,0.65910,0.66276,0.66640,0.67003,0.67364,0.67724,0.68082,0.68439,0.68793},
    {0.5,0.69146,0.69497,0.69847,0.70194,0.70540,0.70884,0.71226,0.71566,0.71904,0.72240},
    {0.6,0.72575,0.72907,0.73237,0.73565,0.73891,0.74215,0.74537,0.74857,0.75175,0.75490},
    {0.7,0.75804,0.76115,0.76424,0.76730,0.77035,0.77337,0.77637,0.77935,0.78230,0.78524},
    {0.8,0.78814,0.79103,0.79389,0.79673,0.79955,0.80234,0.80511,0.80785,0.81057,0.81327},
    {0.9,0.81594,0.81859,0.82121,0.82381,0.82639,0.82894,0.83147,0.83398,0.83646,0.83891},
    {1.0,0.84134,0.84375,0.84614,0.84849,0.85083,0.85314,0.85543,0.85769,0.85993,0.86214},
    {1.1,0.86433,0.86650,0.86864,0.87076,0.87286,0.87493,0.87698,0.87900,0.88100,0.88298},
    {1.2,0.88493,0.88686,0.88877,0.89065,0.89251,0.89435,0.89617,0.89796,0.89973,0.90147},
    {1.3,0.90320,0.90490,0.90658,0.90824,0.90988,0.91149,0.91308,0.91466,0.91621,0.91774},
    {1.4,0.91924,0.92073,0.92220,0.92364,0.92507,0.92647,0.92785,0.92922,0.93056,0.93189},
    {1.5,0.93319,0.93448,0.93574,0.93699,0.93822,0.93943,0.94062,0.94179,0.94295,0.94408},
    {1.6,0.94520,0.94630,0.94738,0.94845,0.94950,0.95053,0.95154,0.95254,0.95352,0.95449},
    {1.7,0.95543,0.95637,0.95728,0.95818,0.95907,0.95994,0.96080,0.96164,0.96246,0.96327},
    {1.8,0.96407,0.96485,0.96562,0.96638,0.96712,0.96784,0.96856,0.96926,0.96995,0.97062},
    {1.9,0.97128,0.97193,0.97257,0.97320,0.97381,0.97441,0.97500,0.97558,0.97615,0.97670},
    {2.0,0.97725,0.97778,0.97831,0.97882,0.97932,0.97982,0.98030,0.98077,0.98124,0.98169},
    {2.1,0.98214,0.98257,0.98300,0.98341,0.98382,0.98422,0.98461,0.98500,0.98537,0.98574},
    {2.2,0.98610,0.98645,0.98679,0.98713,0.98745,0.98778,0.98809,0.98840,0.98870,0.98899},
    {2.3,0.98928,0.98956,0.98983,0.99010,0.99036,0.99061,0.99086,0.99111,0.99134,0.99158},
    {2.4,0.99180,0.99202,0.99224,0.99245,0.99266,0.99286,0.99305,0.99324,0.99343,0.99361},
    {2.5,0.99379,0.99396,0.99413,0.99430,0.99446,0.99461,0.99477,0.99492,0.99506,0.99520},
    {2.6,0.99534,0.99547,0.99560,0.99573,0.99585,0.99598,0.99609,0.99621,0.99632,0.99643},
    {2.7,0.99653,0.99664,0.99674,0.99683,0.99693,0.99702,0.99711,0.99720,0.99728,0.99736},
    {2.8,0.99744,0.99752,0.99760,0.99767,0.99774,0.99781,0.99788,0.99795,0.99801,0.99807},
    {2.9,0.99813,0.99819,0.99825,0.99831,0.99836,0.99841,0.99846,0.99851,0.99856,0.99861},
    {3.0,0.99865,0.99869,0.99874,0.99878,0.99882,0.99886,0.99889,0.99893,0.99896,0.99900},
    {0.0,   0.00,   0.01,   0.02,   0.03,   0.04,   0.05,   0.06,   0.07,   0.08   ,0.09}
};
double percentile_z      = 1.96;               // default for confidence level 0.95
bool   pci               = false;              // default: do not plot proportion confidence interval
std::string   pci_chr    = "";                 // chromosome id  to plot PCI
unsigned long pci_start  = 1;                  // start position to plot PCI
unsigned long pci_end    = 1000;               // end   position to plot PCI
double        pci_cfd    = 0.95;               // confidence level  for  PCI
