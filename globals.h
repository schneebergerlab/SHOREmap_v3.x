// globals.h
#include   <string>
#include <stddef.h>
#include      <map>
#ifndef __globals_h__
#define __globals_h__
#define INF      9999999999
#define MARGIN 0.0000000001

using namespace std;

struct TRIPLE
{
    unsigned long Ci[3];
};

struct WINSAVG // window's average values
{
    unsigned long   wbegin;
    unsigned long   wend;
    double wfreq;
    double wr;
    double wboost;
};

struct LIKELIHOOD // TO REMOVE
{
    double P1est;
    double errEst;
    double min;
};

struct INTERVAL   // TO REMOVE
{
    unsigned long   start;
    unsigned long   size;
    double value; // value
    double lev; // level
};
struct REGION
{
    unsigned long start;
    unsigned long size;
};

struct MARK6
{
    std::string ref;
    std::string mut;
    std::string qua;
    std::string cov;
    std::string ccd;
    std::string aht;
    std::string par;                  // indicates it determined in which parent    
};
////// Command line parameters ///////////////////////////////////////////////////////////////////////////////////////////////////////////
extern unsigned long row_first;       /* the first line to extract from the consensus file               parameters                     */ 
extern unsigned long row_last;        /* the last  line to extract from the consensus file               for                            */
extern unsigned long col_first;       /* the first row  to extract from the consensus file               extracting                     */
extern unsigned long col_last;        /* the last  row  to extract from the consensus file               consensus info                 */
extern unsigned long extract_bg_ref;  /* extracting from quality_reference.txt of a sequenced background mutant                         */
/****************************************************************************************************************************************/
extern char*         mycross;         // outcross or backcross.
extern double        expect;          // expected coverage at a locus.
extern double        interval_min_mean;// minimum mean                     of allele frequency of a set of markers in a mapping interval
extern double        interval_max_mean;// maximum mean                     of allele frequency of a set of markers in a mapping interval
extern double        interval_max_cvar;// maximum coefficient of variation of allele frequency of a set of markers in a mapping interval
extern double        confidence;      // confidence level [default 0.99].
extern std::string   fchrsizes;       // name of a file which contains chromosome id and size in format: <chr chrSize>; to be initialized accoding to input parameters.
extern std::string   out_folder;      // name of a folder where output files can be found; to be initialized accoding to input parameters.
extern std::string   fmarker ;        // name of a file which contains markers in format: <0.project_name 1.chr 2.position 3.ref_base 4.cons_base ...>; note: remianing 4 items <5.quality_score 6.support 7.concordance 8.avg_hits> are not in use.
extern std::string   fmarker_pa;      // name of a file which contains markers of parent a (used for filtering - not necessary)
extern std::string   fmarker_pb;      // name of a file which contains markers of parent b (used for filtering - not necessary)
extern unsigned long chrsizes_max;
extern double        quality_min;     // minimum base quality score (for backcross only)
extern double        quality_max;     // maximum base quality score used in find_marker_pos_2parents_v3.cpp
//marker_format;                      //
extern std::string   fconsensus;      // name of a file which contains alleles in format: <0.chr 1.position 2.cons_base 3.coverage 4.A_support 5.C_support 6.G_support 7.T_support ...>; note: remaining 57 items after the listed are not in use.
//consensus_format;                   //
extern std::string     fconsenref;    // name of a file   which contains alleles	for a reference used for filtering
extern std::string     fpeaks;        // name of a file   which contains peak-positions after calling visualize or backcross function	
extern unsigned long   window_size;          // window size for visualization and peak finding etc [default 50,000].
extern unsigned long   window_step;          // window step for ...[default 10,000].
                                             //
extern unsigned long   peak_window_size;     // window size for peaking finding [default 50,000].
extern unsigned long   peak_window_step;     // window step for ...             [default 10,000].
                                             //
extern unsigned long   filter_min_marker;    // threshold value: if a window has less markers than this minimum, it will be filtered.
extern unsigned long   filter_min_coverage;  // threshold value: if a locus has a coverage less than this minimum, it will be filtered.
extern unsigned long   filter_max_coverage;  // threshold value: if a locus has a coverage more than this maximum, it will be filtered.
extern unsigned long   cluster_avg_coverage; // threshold value: this value will be used for clustering of markers before plotting.
extern unsigned long   cluster_dim_coverage; // threshold value: this value will be used for clustering of markers before plotting as diameter
                                             //
extern unsigned long   outlier_window_size;  // window size to assess local allel frequency used for outlier removal [default 200,000].
extern double          outlier_pvalue;       // p_value for outlier removal.
                                             //
extern double          misphenotyped;        // degree of putatively mis-scored plants [default 0.00]. ???
                                             //
extern std::string     reg_chromosome;       // name/id of the chromosome to zoom. [!!!caution: SHOREmap.pl cannot recognize string names/ids].
extern unsigned long   reg_begin;    // to zoom from this point.
extern unsigned long   reg_end;      // to zoom   to this point.
extern double          reg_freq_min; // minimum frequency of zoom region.
extern double          reg_freq_max; // maximum frequency of zoom region.
                                     //
extern std::string      freferror;   // name of a file which contains reference errors.
extern unsigned long    background2; // the mutation is in the second parent if it equals 1 [default 0].        
extern unsigned long    verbose;     // to be talktive or not during process.
                                     //
extern unsigned long    runid;       // running id of a specific task.
extern unsigned long    r_max;       // distance of a specific mutation from a peak.
extern bool             plot_r;      // to plot frequency calculation 'r' instead of 'boost'.
extern double           boost_max;   // maximum boost value.
extern bool             plot_boost;  // to plot 'boost'.
extern bool             plot_marker; // to add single markers in the plot.
extern unsigned long    plot_bc;;    // switch on or not plotting of markers from backcross process
extern bool             plot_window; // plot window-averaged allele frequency of markers
extern bool             plot_scale;  // plot the chromosomes scaled to the one with the maximum length
extern unsigned long    plot_record; // record statsitcal info that has been plotted in pdf
                                     //
extern unsigned long    marker_score;// minimum score cutoff for filtering markers [default 25].
//extern unsigned long  marker_read; // minimum read support [default 0].
extern double           marker_freq; // minimum concordance: ratio of reads supporting a predicted feature to total coverage [default 20].      
extern double           marker_hit;  // maximum average number of hit about a marker  
extern unsigned long    pmarker_score;   // minimum score cutoff of parent markers for filtering markers [default 0 - no par filtering].                    
extern unsigned long    pmarker_min_cov; // minimum coverage cutoff of parent markers for filtering markers [default 0 - no par filtering].                  
extern unsigned long    pmarker_max_cov; // maximum coverage cutoff of parent markers for filtering markers [default INF - no par filtering]. 	              
extern double           pmarker_min_freq;// minimum frequency of parent markers for filtering markers [default 0.0 - no par filtering]   
extern std::string      pmarker_ab_ratio;// ratio of markers of parents [default 1,1,1,1 - the cutoffs for parent-markers are the same]               
extern std::string      frun_file;   // to indicate which marker file is in use.
extern unsigned long    indel_size;  // size of indels to extract from a vcf file
                                     //
/* filtering according to ref info   */
extern unsigned long    bg_ref_filter;// trun on this filter or not
extern std::string      bg_ref_base_file;
extern std::string      bg_ref_base_file_pb;
extern unsigned long    bg_ref_cov;
extern unsigned long    bg_ref_cov_max;
extern double           bg_ref_freq;
extern unsigned long    bg_ref_score;
extern unsigned long    fg_N_cov;    // <=> marker_N_cov  
extern unsigned long    fg_INDEL_cov;// <=> marker_INDEL_cov
                                     //
extern std::string fbackground_file; // name of files which contain backgroud mutations (!!!caution: comma saparated if more than one).
extern unsigned long    bg_score;    // minimum score cutoff for filtering background markers [default 0].
extern unsigned long    bg_read;     // minimum read support for background markers [default 0].
extern double           bg_freq;     // minimum concordance for background markers [default 20].

//////////ploting related vari////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern unsigned long    clusterK;      // levels of scores for ranking markers in visualization [default 40]
extern unsigned long    summary;       // 
extern unsigned long    only_EMS;      // only ploting EMS-induced mutations 
extern unsigned long    other_mutant;  // other mutagen is not used for screening
extern unsigned long    filter_plot;   // to plot after filtering
extern std::string      R_input;       // name of input when calling R to plot [caution: this is removed in this version using C].
extern std::string      R_out;         // folder of output when ...
extern int              vis_annotation;// visualization of annotation along chromosomes

////// Additional global variables   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

extern map<std::string, unsigned long>  CHR2SIZE; // where chromosome and its size is recorded: CHR2SIZE["chr"]     = (unsigned long)size.
extern map<std::string, bool>  REFERROR;          // this marker is an error:                   REFERROR["chr.pos"] = 1.
extern map<std::string, std::string>  ALLELE1;    // where the mutated allele is recorded:      ALLELE1["chr.pos"]  = 'A/C/G/T'. - af to calculate
extern map<std::string, std::string>  ALLELE2;    // where the ref allele is recorded:          ALLELE2["chr.pos"]  = 'A/C/G/T'.
extern map<std::string, std::string>  QUALITY1;   // where quality of base call of the mutated allele is recorded: QUALITY1["chr.pos"]
extern map<std::string, std::string>  CMD;        // where the command line parameters are recorded.  
extern std::string  strcatCMD;                    // a string concatenation of the cmds
extern std::string  pwd;                          // path of working folder.
extern map<std::string, std::string> bgMARKER;    // information about bg-markers which have been used to filtering/retaining foreground markers
extern map<std::string, std::string> bgREF;       // information about bg reference bases which have been used to keep a fg-marker
extern map<std::string, std::string> bgREFB4A;    // information about bg reference bases according to parent B which have been used to keep a fg-marker in parent A - used in outcross clustering
extern map<std::string, std::string> bgREFA4B;    // information about bg reference bases according to parent A which have been used to keep a fg-marker in parent B - used in outcross clustering
extern bool keep_common;                          // whether keep common snps of parental lines in marker list - create - default false

////// Collect marker counts and sliding window values /////////////////////////////////////////////////////////////////////////////////////////////////////////

extern map<std::string, map<std::string, unsigned long> > CHR2POS2ALLELE1_COUNT; // where the count for allele 1 is recorded: CHR2POS2ALLELE1_COUNT["chr"][pos] = allele1_count.
extern map<std::string, map<std::string, unsigned long> > CHR2POS2ALLELE2_COUNT; // where the count for allele 2 is recorded: CHR2POS2ALLELE2_COUNT["chr"][pos] = alelle2_count.
extern map<std::string, map<std::string, unsigned long> > CHR2POS2ERROR_COUNT;   // where the count for error    is recorded: CHR2POS2ERROR_COUNT["chr"][pos]   = error_count.
extern map<std::string, map<unsigned long, TRIPLE>  >     CHR2POS2_ale1_ale2_err_COUNT; /* shq: postion records like allele 1, allele 2, error - 2013-03-10 19:31   */

/* file for markers annotation */
extern std::string snp_file;
extern std::string del_file;
extern std::string ins_file;
extern std::string gff_file;
extern std::string genome_file;
/* info about a single SNP */
struct mySNP
{
    std::string   ecotype;
    std::string   chromosome;
    unsigned long position;
    std::string   stype;
    std::string   gene_id;
    unsigned long gene_pos;  // the position of CDS relative to the start/end of a gene
    unsigned long cds_pos;   // the relative position of CDSs if counted from the start of a gene, e.g., 1, 2, 3, cds_length...
    unsigned long codon_pos; // in {1,2,3}: position of a SNP on a codon
    unsigned long ns_change;
    unsigned long new_stop;
    unsigned long lost_stop;
    unsigned long splicechange;
    std::string   ref_base;
    std::string   new_base;
    std::string   ref_aa;
    std::string   new_aa;
    //map<std::string, std::string>domain_change;
    unsigned long support;
    double        concordance;
    double        quality;
    unsigned long peak_distance;
    double        marker_ratio;
    std::string   source;
};
extern multimap<unsigned long, mySNP> SNPlist; // SNP list to be annotated with gene info (gff file)

/* info about an insertion or deletion */
struct myINDEL
{
    std::string   ecotype;
    std::string   chromosome;
    unsigned long begin;
    unsigned long end;
    std::string   seq;       // seq of indel
    std::string   stype;
    std::string   gene_id;
    unsigned long gene_pos;  // the position of CDS relative to the start/end of a gene
    unsigned long cds_pos;   // the relative position of CDSs if counted from the start of a gene, e.g., 1, 2, 3, cds_length...
    unsigned long codon_pos; // in {1,2,3}: position of a SNP on a codon
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
extern multimap<unsigned long, myINDEL> INSlist;  // INSERTION list to be annotated with gene info (gff file)
extern multimap<unsigned long, myINDEL> DELlist;  // DELETION  list to be annotated with gene info (gff file)

/* info about gene annotation */
struct QUARTET                                    // it is not quartet anymore as frame info is added -- 2016-MAR-03.
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
extern map<std::string, QUARTET> gene_ann;
extern map<std::string, map<unsigned long, QUARTET> > coding_ann;
extern multimap<std::string, multimap<std::string, RANGE> > seq_type; // <gene_name, <"CDS/gene/...", RANGE> >
extern int first_codon_frame;                     // frame value for the first codon: can be '0','1','2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
extern int last_codon_frame;                      // frame value for the last codon: can be '0','1','2'.

/* info about GeneSNPlist */
struct GeneSNPlist
{
    std::string   gene_id;                        //1
    std::string   isoform;                        // 
    std::string   chromosome;                     //
    unsigned long start;                          //
    unsigned long end;                            //5                
    unsigned long cds_length;                     //           
    unsigned long protein_length;                 //
    std::string   orientation;                    //
    multimap<unsigned long, mySNP> SNP_lists;     //
    map<unsigned long, unsigned long> coding_SNP; //10 first=virtual(1,2,3,...,cds_length), second=real(along genome)
    map<unsigned long, unsigned long> NS_changes; //
    map<unsigned long, unsigned long> AA_changes; //
    map<unsigned long, unsigned long> CDSexon;    //
    map<std::string, std::string> protein;        //
    std::string  ref_gene;                        //15
    std::string  ref_coding;                      //
    std::string  eco_coding;                      //17
};

/* log with hints for solving possible errors     */
extern std::string errLog;  // TODO 2013-04-19 09:56
/* z values for calculating confidence interval   */
extern double z_table[32][11];
extern double percentile_z;
extern bool   pci;                                // default: do not plot proportion confidence interval
extern std::string   pci_chr;
extern unsigned long pci_start;
extern unsigned long pci_end;
extern double        pci_cfd;
#endif
