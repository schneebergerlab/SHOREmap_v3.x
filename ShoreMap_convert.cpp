/* This function converts the resequencing results in vcf format to SHOREmap-acceptable format    
   if AF>0.2, will be recorded as quality_reference.txt
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include  <fstream>
#include <iostream>
#include <vector>

// file, dir
#include    <dirent.h>
#include  <sys/stat.h>
#include <sys/types.h>
#include    <unistd.h>
#include      <time.h>

#include "globals.h"
#include "split_string.h"
#include "precheck_opt.h"
#include "print_error_exit.h"
#include "convert_shore_pileup.h"

bool vcf2shoremap(char* VCFmarkerfile, unsigned long indelsize, int runid, char* out_path);
std::string get_indelseq_bycmp_refalt(std::string sref, std::string salt);
int indel2snp(std::string iref, std::string ialt, std::string* sref, std::string* salt, int* offset);

bool record_snp = true;
bool record_con = true;
bool record_ref = true;
bool cpileup    = false;

bool ShoreMap_convert(int argc, char* argv[])
{
    if(argc < 3)
    {
        printf("\nThis function converts info: in VCF4.1 format into one accepted by SHOREmap\n");
        printf("                          or in shore pileup into shore consensus_summary.   \n\n");
        printf("Usage: SHOREmap convert [options]                                            \n");
        printf("#Mandatory:                                                                  \n");
        printf("\t--marker     STRING     VCF file of variations or shore pileup file from qVar\n");
        printf("\t--folder     STRING     output folder.                                     \n");
        printf("#Optional:                                                                   \n");
        printf("\t--indel-size INT        Size of indels [0] (if 0, cancel output indels).   \n");
        printf("\t--min-AF     DOUBLE     vcf*: minimum AF to record quality reference [0.2] \n");
        printf("\t-no-v        BOOL       vcf*: do not record snp variants [FALSE]           \n");
        printf("\t-no-c        BOOL       vcf*: do not record all consensus bases [FALSE]    \n");
        printf("\t-no-r        BOOL       vcf*: do not record quality reference bases [FALSE]\n");
        printf("\t-p           BOOL       when it is on, converting shore pileup file [FALSE]\n");
        printf("\t-runid       INT        ID of run [1]                                      \n");
        printf("vcf*: option only works for coverting files in VCF    \n");
        printf("\nExample1: SHOREmap convert --marker samtools.vcf --folder testConvert -runid 4;");
        printf("\nExample2: SHOREmap convert --marker pileup.txt -p --folder testConvert -runid 5.\n\n");
        exit(1);
    }
    
    double startT=clock();
    /* set cmd line inputs .......................................................................*/
    runid        = 1;
    indel_size   = 0;
    reg_freq_min = 0.2;
    int ic = 2;
    while (ic < argc) 
    {
        strcatCMD += " ";
        strcatCMD += (string)argv[ic];
        if(!strcmp(argv[ic],"--marker"))     // option 1 to variable: fmarker.......................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, false);
            fmarker        += argv[ic];
            FILE* fp_marker = fopen((char*)fmarker.c_str(), "r"); 
            if(fp_marker == NULL)
            {
                printf("marker file \'%s\' does NOT exist. Exited.\n", fmarker.c_str());
                exit(1);
            }
            printf("\tFile to be converted:\t\t%s\n", fmarker.c_str());
           // read_marker();
            fclose(fp_marker);
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
                 /* if !mkdir() is TRUE: a new directory has to be created.
                    if out_folder is "outfolder/" instead of "/your/output/path/outfolder", 
                    an "outfolder/" will be created under the current working directory */
                    printf("\tFolder\t\t\t\t\t%s created.\n", out_folder.c_str());
                }
                else 
                {
                    printf("ERROR: cannot create output path. Exited.\n");
                    exit(1);
                }
            }
            closedir(dir_out_folder);
        }
        else if(!strcmp(argv[ic],"--indel-size"))  // option 3 to variable: indel_size..............
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            indel_size = atol(argv[ic]);
            if(indel_size<0) 
                { printf("ERROR: arg of %s must be >= 0. Exited.\n",
                       argv[ic-1]); exit(1);}
            printf("\tindel_size set as:\t\t%ld\n", indel_size);
        }
        else if(!strcmp(argv[ic],"--min-AF"))     // option 3.5 to variable: reg_freq_min..........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            reg_freq_min     = atof(argv[ic]);
            if(reg_freq_min<0.0 || reg_freq_min>1.0)
                { printf("ERROR: arg of %s must be [0.0, 1.0]. Exited.\n", argv[ic-1]); exit(1);}
            printf("\treg_freq_min set as:\t\t\t%.2f.\n", reg_freq_min);
        }
        else if(!strcmp(argv[ic],"-runid"))        // option 4 to variable: runid...................
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], true, true);
            runid = atoi(argv[ic]);
            printf("\trunid set as:\t\t\t\t%ld.\n", runid);
        }
        else if(!strcmp(argv[ic],"-no-v"))        // option 5 to variable: record_snp........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            record_snp         = false;
            printf("\trecord_snp set as:\t\t\t%d (off snp-output)\n", record_snp);
        }
        else if(!strcmp(argv[ic],"-no-c"))        // option 6 to variable: record_con........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            record_con         = false;
            printf("\trecord_con set as:\t\t\t%d (off consensus-output)\n", record_con);
        }
        else if(!strcmp(argv[ic],"-no-r"))        // option 7 to variable: record_ref........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            record_ref         = false;
            printf("\trecord_ref set as:\t\t\t%d (off quality-reference-output)\n", record_ref);
        }
        else if(!strcmp(argv[ic],"-p"))        // option 8 to variable: cpileup........
        {
            ic = precheck_opt(argc, argv, ic, argv[ic], false, false);
            cpileup         = true;
            printf("\tcpileup set as:\t\t\t%d (converting shore pileup file)\n", cpileup);
        } 
        else // other options not necessary
        { 
            printf("\tWarning: \"%s\" is NOT a valid option.\n", argv[ic]);
        }
        ic ++;
    }
    if(CMD.find("--marker") == CMD.end())
    {
        print_error_exit((char*)"--marker", false);
    }
    else 
    if(CMD.find("--folder") == CMD.end())
    {
        print_error_exit((char*)"--folder", false);
    }
    
    if(!cpileup && !record_snp && !record_con && !record_ref)
    {
        printf("\tERROR: the converted information to record must be provided. \n");
        printf("\tHINT: pls remove at least one option of -no-v, -no-c, and -no-r from your cmd.\n");
        printf("\tExited!\n\n");
        exit(1);
    }
    
    if(cpileup && convert_shore_pileup(runid, (char*)fmarker.c_str()) > 0)
    {
        /* convert shore pileup to consensus_summary.txt format */
        printf("ERROR: shore pipleup not correctly converted. Check the files.\n");
        exit(1);
    }
    else if(!cpileup && !vcf2shoremap((char*)fmarker.c_str(), indel_size, runid, (char*)out_folder.c_str()))
    {
        /* convert VCF to SHOREmap acceptable format */
        printf("ERROR: VCF not correctly converted. Check the files.\n");
        exit(1);
    }
    
    double finishT = clock();
    printf("\nTime consumed %.4f seconds.\n", (finishT-startT)/1000000);
   
    return true;
}

/* function to convert vcf to SHOREmap acceptable format */
bool vcf2shoremap(char* VCFmarkerfile, unsigned long indelsize, int runid, char* out_path)
{
    /* input file */
    std::ifstream fpin (VCFmarkerfile);
    if(!fpin.is_open())
    {
        printf("Marker file \'%s\' does NOT exist. Exited.\n", VCFmarkerfile);
        exit(1);
    }
    printf("\tReading marker read from vcf file:\t%s...\n", VCFmarkerfile);
    
    /* output file: snp ..........................................................................*/
    char convertedmkrfile[1024];
    sprintf(convertedmkrfile, "%s%d_converted_variant.txt", out_path, runid);
    FILE* fpout;
    if(record_snp)
    {
        fpout = fopen(convertedmkrfile, "w");
        if(fpout == NULL)
        {
            printf("Cannot open file to write converted snp info.\n");
            exit(1);
        }
    }
    /* output file: consensus ....................................................................*/
    char convertedconsenfile[1024];
    sprintf(convertedconsenfile, "%s%d_converted_consen.txt", out_path, runid);
    FILE* fpoutCON;
    if(record_con)
    {
        fpoutCON = fopen(convertedconsenfile, "w");
        if(fpoutCON == NULL)
        {
            printf("Cannot open file to write converted consensus info.\n");
            exit(1);
        }
    }
    /* output file: quality_reference ............................................................*/
    char convertedreferencefile[1024];
    sprintf(convertedreferencefile, "%s%d_converted_reference.txt", out_path, runid);
    FILE* fpoutREF;
    if(record_ref)
    {
        fpoutREF = fopen(convertedreferencefile, "w");
        if(fpoutREF == NULL)
        {
            printf("Cannot open file to write converted quality reference info.\n");
            exit(1);
        }
    }
    
    // TODO: get indel size from cmd line; 0 means does not output indels
    /* output file: insertion */
    char convertedinsfile[1024];
    FILE* fpoutINS;
    if(indelsize > 0)
    {
        sprintf(convertedinsfile, "%s%d_converted_insertion_info.txt", out_path, runid);
        fpoutINS = fopen(convertedinsfile, "w");
        if(fpoutINS == NULL)
        {
            printf("Cannot open file to write converted insertion info.\n");
            exit(1);
        }
    }
    /* output file: deletion */
    char converteddesfile[1024];
    FILE* fpoutDEL;
    if(indelsize > 0)
    {
        sprintf(converteddesfile, "%s%d_converted_deletion_info.txt", out_path, runid);
        fpoutDEL = fopen(converteddesfile, "w");
        if(fpoutDEL == NULL)
        {
            printf("Cannot open file to write converted deletion info.\n");
            exit(1);
        }
    }
    
    unsigned long num_snp = 0;
    unsigned long num_ins = 0;
    unsigned long num_del = 0;
    unsigned long num_skp = 0;
    unsigned long num_con = 0;
    unsigned long num_snp_line = 0;
    const char* token;
    std::string line;
    while(fpin.good())
    {
        line.clear();
        getline(fpin, line);
        if (line.length() == 0 || line.find("#") != std::string::npos) {continue;} // # header
        
        std::vector<std::string> splittedLine = split_string(line, '\t');
        if(splittedLine.size() < 8)
        {
            /* the first 8 columns are mandatory in vcf4.1 */
            printf("\tInsufficient info @ %s.\n", line.c_str());
            printf("\tHint: at least 8 columns required.\n");
            printf("\tSkipped.\n");
            continue;
        }
        // splittedLine[0]: CHROM
        // splittedLine[1]: POS
        // splittedLine[2]: ID
        // splittedLine[3]: REF
        // splittedLine[4]: ALT
        // splittedLine[5]: QUAL
        // splittedLine[6]: FILTER
        // splittedLine[7]: INFO
        // splittedLine[8]: FORMAT
        // splittedLine[9]: VALUE

        int necessaryCheck    = 0;
        double        iAF     = 0;  // alt allele frequency
        std::string   aDP     = "0";// overall raw coverage
        unsigned long iDP     = 0;  // read coverage of ref+mut
        unsigned long A_cov   = 0;  
        unsigned long C_cov   = 0;
        unsigned long G_cov   = 0;
        unsigned long T_cov   = 0;
        unsigned long D_cov   = 0;  // caution: D_cov and N_cov is not extractable from a vcf file.
        unsigned long N_cov   = 0;  // caution:
        unsigned long ref_cov = 0;
        unsigned long alt_cov = 0;
        std::vector<std::string> splittedInfo = split_string(splittedLine[7], ';');
        // caution: avoid values with symbol '*' - TODO
        bool givenDP  = false;
        bool givenDP4 = false;
        bool givenAF  = false;
        for (int si = 0; si < splittedInfo.size(); si++) 
        {
            std::vector<std::string> miniInfo = split_string(splittedInfo[si], '=');
            if(miniInfo[0] == "DP")
            {
                aDP     = miniInfo[1];
                //ref_cov = atol(aDP.c_str()); // to remove
                //alt_cov = 0;                 // to remove
                necessaryCheck ++;    
                givenDP = true;            
            }
            else 
            if(miniInfo[0] == "DP4")
            {
                // like DP4=0,43,0,38 (first two: ref-coverage, last-two: alt coverage)
                std::vector<std::string> DP4Info = split_string(miniInfo[1], ',');
                ref_cov  = atol(DP4Info[0].c_str())+atol(DP4Info[1].c_str());
                alt_cov  = atol(DP4Info[2].c_str())+atol(DP4Info[3].c_str());
                iDP      = alt_cov + ref_cov;
                necessaryCheck ++;
                givenDP4 = true;
            }
        }
        if(givenDP && givenDP4)
        {
            /* caution: raw AF with error - 2014-05-05 */
            if(atol(aDP.c_str()) <= 0)                                  // suggested by Emily Wynn - University of Nebraska–Lincoln
            iAF     = (double)alt_cov/(double)atol(aDP.c_str());        // caution: in case no coverage, divided by 0 - should never happen!
            else
            iAF     = (double)alt_cov/(double)iDP;                      // this always happens!
        }
        else
        if(splittedLine[8].compare("GT:AD:DP:GQ:PL")==0)
        {
            // splittedLine[8]: GT:AD:DP:GQ:PL ==> splittedLine[9]: 0/1:108,19:127:99:463,0,4548
            vector<string> sampleinfo = split_string(splittedLine[9], ':');
            vector<string> ADinfo     = split_string(sampleinfo[1], ',');
            if(ADinfo.size()<2) 
            {
                printf("\tWARNING: Unexpected INFO @line (GT:AD:DP:GQ:PL, AD requires 2 values) ): %s\n", (char*)line.c_str());
            }
            ref_cov = atol(ADinfo[0].c_str());
            alt_cov = atol(ADinfo[1].c_str());
            if(ADinfo.size()==2)
            {
                iDP     = ref_cov + alt_cov;
            }
            else
            {
                iDP     = atol(sampleinfo[2].c_str());
            }
            iAF         = (double)alt_cov/(double)iDP;
            // for vcf4.2
            aDP         = sampleinfo[2];
        }
        else
        {
            printf("\tWARNING: Unexpected INFO @line (either (DP, DP4) or (GT:AD:DP:GQ:PL) required): %s\n", (char*)line.c_str());
            continue;
        }
        
        // check if it is an indel: based on alt allele
        bool var_is_indel = false;
        vector<string> altinfo = split_string(splittedLine[4], ',');
        // caution: only consider the major allele
        if(altinfo[0].size() != splittedLine[3].size()) 
        {
            var_is_indel = true;
        }
        
        
        if(indelsize>0 && (line.find("INDEL")!=std::string::npos || var_is_indel==true))
        {
             // CAUTION: "IS" info is not as general as "DP4", but more accurate.
             // length of indels considered - TODO?
             if(splittedLine[4].find(".")!=std::string::npos && splittedLine[4].size()==1)
             {
                /* if alleleInfo==".", which means no alt allele; how do we handle such INDELS?   */
                ////printf("WARNING: unexpected INDEL @line: %s.\n", line.c_str());     // caution!!
                continue;
             }
             /* Caution!!! only consider the major/first alt allele ..............................*/
             std::vector<std::string> indelseq4Info = split_string(splittedLine[4], ',');
             if(indelseq4Info.size() > 1) 
             {
                 printf("Warning: multiple alt-alleles at line %s (only major allele extracted).\n", line.c_str());
             }
             if(splittedLine[3].size() < indelseq4Info[0].size()) // ref shorter => insertion in mut
             {
                 if(indelseq4Info[0].size() - splittedLine[3].size() > indelsize) 
                 {
                     num_skp ++;
                     continue; 
                 }
                 /* CAUTION: only major insertion allele */
                 std::string seq_ins = get_indelseq_bycmp_refalt(splittedLine[3], indelseq4Info[0]);
                 /* write insertion info */
                 if(num_ins > 0) fprintf(fpoutINS, "\n");
                 fprintf(fpoutINS, "proj_%d\t", runid);           // 0.proj name
                 fprintf(fpoutINS, "%s\t%s\t%ld\t%ld\t%s\t%s\t%ld\t%.4f\t%.4f", 
                         (char*)splittedLine[0].c_str(),          // 1.chr
                         (char*)splittedLine[1].c_str(),          // 2.start
                         atol(splittedLine[1].c_str())+1,         // 3.end
                         seq_ins.length(),                        // 4.length
                         (char*)seq_ins.c_str(),                  // 5.ins_seq
                         (char*)"UNKWN",                          // 6.read_type
                         //iDP,                                   // not used since 2013-06-03
                         alt_cov,                                 // 7.support - caution - only mut
                         iAF,                                     // 8.concordance
                         1.0000);                                 // 9.avg_hits
                 num_ins ++;
            }
            else if(splittedLine[3].size() > indelseq4Info[0].size()) // ref longer => deletion in mut
            {
                if(splittedLine[3].size() - indelseq4Info[0].size() > indelsize)
                {
                     num_skp ++;
                     continue; 
                }
                /* CAUTION: only major deletion allele */
                std::string seq_ins = get_indelseq_bycmp_refalt(splittedLine[3], indelseq4Info[0]);
                /* write deletion info */
                if(num_del > 0) fprintf(fpoutDEL, "\n");
                fprintf(fpoutDEL, "proj_%d\t", runid);           // 0.proj name
                fprintf(fpoutDEL, "%s\t%ld\t%ld\t%ld\t%s\t%s\t%ld\t%.4f\t%.4f", 
                        (char*)splittedLine[0].c_str(),          // 1.chr
                        atol(splittedLine[1].c_str()) + 1,       // 2.start
                        atol(splittedLine[1].c_str())+seq_ins.length(), // 3.end
                        seq_ins.length(),                        // 4.length
                        (char*)seq_ins.c_str(),                  // 5.ins_seq
                        (char*)"UNKWN",                          // 6.read_type
                        //iDP,                                   // not used since 2013-06-03
                        alt_cov,                                 // 7.support - caution - only mut
                        iAF,                                     // 8.concordance
                        1.0000);                                 // 9.avg_hits
                num_del ++;
            }
            else if(splittedLine[3].size() == indelseq4Info[0].size())
            {
                /* the major allele indicates SNP change */
                printf("\tINDEL to SNP line (discarded): %s\n", line.c_str());
                
                /* TODO 2014-01-12: handle this kind of INDEL case: ref=ACG, alt=ATG,AC */
                std::string sref("");
                std::string salt("");
                int offset = 0;
                indel2snp(splittedLine[3], indelseq4Info[0], &sref, &salt, &offset);
            }
        }
        else // SNPs
        if(line.find("INDEL") == std::string::npos && var_is_indel==false)
        {
            num_snp_line ++;
            /* ref allele coverage */
            if(splittedLine[3]=="A") A_cov = ref_cov; else 
            if(splittedLine[3]=="C") C_cov = ref_cov; else 
            if(splittedLine[3]=="G") G_cov = ref_cov; else 
            if(splittedLine[3]=="T") T_cov = ref_cov; else
            {
                ////printf("WARNING: SNP - unexpected ref base @line: %s.\n", line.c_str());
                ////continue;
            }
            /* alt allele coverage */
            std::vector<std::string> alleleInfo;
            if(splittedLine[4].find(".")!=std::string::npos && splittedLine[4].size()==1)
            {
                /* if alleleInfo==".", which means no alt allele; we assume alt=ref               */
                alleleInfo.push_back(splittedLine[3]);
                if(alleleInfo[0]=="A")  A_cov = ref_cov; else 
                if(alleleInfo[0]=="C")  C_cov = ref_cov; else 
                if(alleleInfo[0]=="G")  G_cov = ref_cov; else 
                if(alleleInfo[0]=="T")  T_cov = ref_cov; else
                {
                    ////printf("WARNING: SNP - unexpected alt base @line: %s.\n", line.c_str());
                }
            }
            else
            {
                alleleInfo       = split_string(splittedLine[4], ',');
                if(alleleInfo.size() > 0)
                {
                    if(alleleInfo[0]=="A")  A_cov = alt_cov; else 
                    if(alleleInfo[0]=="C")  C_cov = alt_cov; else 
                    if(alleleInfo[0]=="G")  G_cov = alt_cov; else 
                    if(alleleInfo[0]=="T")  T_cov = alt_cov; else
                    {
                        ////printf("WARNING: SNP - unexpected alt base @line: %s.\n", line.c_str());
                    }
                }
            }
            /* record the count covering the alternative alleles as errors - recording on 2nd allele
                e.g., ref=A, ale1=C, ale2=G, ale3=T, then A=ref_cov, C=alt_cov, alt ales G+T=aDP-iDP
               else record the unknown counts as number of N calls.
            */
            if(alleleInfo.size() > 1)
            {
                if(alleleInfo[1]=="A") A_cov += atol(aDP.c_str()) - iDP; else 
   	        if(alleleInfo[1]=="C") C_cov += atol(aDP.c_str()) - iDP; else 
        	if(alleleInfo[1]=="G") G_cov += atol(aDP.c_str()) - iDP; else 
        	if(alleleInfo[1]=="T") T_cov += atol(aDP.c_str()) - iDP; else
        	{
        	    ////printf("WARNING: unexpected 2nd alt base @line: %s.\n", line.c_str());
        	}
            }
            //else
            //{
                N_cov = atol(aDP.c_str()) - (A_cov + C_cov + G_cov + T_cov + D_cov);
            //}
            
            
            /* write snp info .................................................................*/
            if(record_snp && num_snp > 0)
            {
                 /* new line */
                 if(alleleInfo[0] != splittedLine[3]) fprintf(fpout, "\n");
            }
            if(record_con && num_con > 0)
            {
                 fprintf(fpoutCON, "\n");
            }
            /* write snp file info */
            if(record_snp && alleleInfo[0] != splittedLine[3])
            {
                if(alleleInfo.size()>=2)
                {
                    printf("Warning: multiple alt-alleles at line %s (only major allele extracted).\n", line.c_str());
                }
                for(int ialt = 0; ialt < 1; ialt ++) // in case multiple alleles
                {
                    if(ialt > 0) fprintf(fpout, "\n");
                    fprintf(fpout, "proj_%d_alt%d\t", runid, ialt);  // 0.proj name
                    fprintf(fpout, "%s\t%s\t%s\t%s\t%s\t",
                     (char*)splittedLine[0].c_str(),
                     (char*)splittedLine[1].c_str(),
                     (char*)splittedLine[3].c_str(),
                     (char*)alleleInfo[ialt].c_str(),               // caution: all alleles recorded
                     (char*)splittedLine[5].c_str());               // 1.chr, 2.pos, 3.ref, 4.alt, 5.qual
                    //fprintf(fpout, "%ld\t%.4f\t1", iDP, iAF);     // 6.ref+alt cov, 7.alt-allele-freq, 8.1
                    fprintf(fpout, "%ld\t%.4f\t1", alt_cov, iAF);   // 6.alt cov,     7.alt-allele-freq, 8.1
                    num_snp ++;
                }
            }
            
            /* wirte quality reference base info .................................................*/
            double ref_AF = (double)ref_cov/(double)atol(aDP.c_str());//caution with error!-2014-05-05
            if(record_ref && ref_AF > reg_freq_min)
            {
                fprintf(fpoutREF, "proj_%d\t", runid);        // 0.proj name  
                fprintf(fpoutREF, "%s\t%s\t%s\t%s\t%s\t",
                 (char*)splittedLine[0].c_str(),
                 (char*)splittedLine[1].c_str(),
                 (char*)splittedLine[3].c_str(),
                 (char*)splittedLine[3].c_str(),              //
                 (char*)splittedLine[5].c_str());             // 1.chr, 2.pos, 3.ref, 4.alt,  5.qual
                fprintf(fpoutREF, "%ld\t%.4f\t1\n", ref_cov, ref_AF);// 6.alt cov,
                                                                     // 7.alt-allele-freq,8.1
            }
            /* write consensus info   */
            if(record_con)
            {
                fprintf(fpoutCON, "%s\t%s\t%s\t%s\t",
                        (char*)splittedLine[0].c_str(),
                        (char*)splittedLine[1].c_str(),
                        (char*)alleleInfo[0].c_str(),               // caution: all alleles recorded
                        (char*)            aDP.c_str());            // 0.chr, 1.pos, 2.alt, 3.raw-cov
                fprintf(fpoutCON, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t1", 
                         A_cov, C_cov, G_cov, T_cov, D_cov, N_cov); // 4.A, 5.C, 6.G, 7.T, 8.Del, 9.N, 
                                                                    // 10. avg_hit=1(fake default)
            }
            num_con ++;
            // OTHER INFO?
        }
    }
    fpin.close();
    if(record_snp)
    fclose(fpout);
    if(record_con)
    fclose(fpoutCON);
    if(record_ref)
    fclose(fpoutREF);
    // TODO: write it into log
    printf("\n\t%ld of (raw_num=%ld) snps, %ld consensus info ", num_snp, num_snp_line, num_con);
    if(indelsize > 0)
    {
        fclose(fpoutINS);
        fclose(fpoutDEL);
        printf(", %ld insertions, %ld deletions (size < %ld; ", num_ins, num_del, indelsize);
        printf("while %ld indels of larger sizes skipped) ", num_skp);
    }
    printf("have been read.\n");
    printf("\t*Note1: there can be positions with reference base as N in the recorded files.\n");
    printf("\t*Note2: unknown counts of reads are recorded as N calls (in last-second column of consensus file). \n");
    return true;
}

/* this function compares ref and major ale sequences given in vcf file to get the pure INDEL seq */
std::string get_indelseq_bycmp_refalt(std::string sref, std::string salt)
{
    /* cases:
     	1: ref=TNN                    alt=TANN                        - insert after T
     	2: ref=GTTT                   alt=GTTTT                       - insert after G
     	3: ref=TCG                    alt=TCAG                        - insert after TC??? - not handled yet!!
     	4: ref=TTATATATATATATATATATAT alt=TTATATATATATATATATAT        - delete after 1st T
     	5: ref=ATC                    alt=A                           - delete after A
     	6: ref=TCG                    alt=TG                          - delete after T
     	so deletion always starts after the first reference base;
     	
    */
    /* NOTE:
       According to vcf4.1 (IGV visualize bam), indels always happen from 2nd position of long seq, right? CHECK!
       That is, first bases of ref and alt seqs are always the same.
       can be removed after future verification, if exceptions happen!!!
    */
    if(sref[0] != salt[0])
    {
        printf("WARNING: skipped indel in vcf: ref=%s\tmut=%s? \n", 
              (char*)sref.c_str(), (char*)salt.c_str());
        //exit(1); // caution here -- 2014-JAN-10 --------> check the case of 'ref=ATC	mut=.'
    }
    std::string slong("");
    std::string sshort("");
    std::string sindel("");
    if(sref.length() > salt.length())
    {
        slong  = sref;
        sshort = salt;
    }
    else 
    if (sref.length() < salt.length())
    {
        slong  = salt;
        sshort = sref;
    }
    else
    {
        printf("Same length of ref and alt sequence in vcf file, but annotated as INDEL, case not supported; ");
        printf("missed in this conversion. Check!\n");
        printf("seq: %s\t%s. Exited (in ShoreMap_convert(...)).\n", 
              (char*)sref.c_str(), (char*)salt.c_str());
        exit(1);
    }
    unsigned long iextr     = 1; // first reference base; then the indel starts??
    // find out the position where ref seq and alt seq start to differ
    //while(iextr < sshort.length())
    //{
    //    if(slong[iextr] != sshort[iextr]) break;
    //    iextr ++;
    //}
    
    unsigned long indel_len = slong.length() - sshort.length(); // length of deleted or inserted

    unsigned long extr_len = 0;
    while(extr_len < indel_len)
    {
        sindel   += slong.substr(iextr, 1);
        iextr    ++;
        extr_len ++;
    }
    return sindel;
}

int indel2snp(std::string iref, std::string ialt, std::string* sref, std::string* salt, int* offset)
{
    int i;
    for(i=0; i<iref.size(); i++)
    {
        if(iref[i] != ialt[i])
        {
            break;
        }
    }
    (*sref)  += iref[i];
    (*salt)  += ialt[i];
    
    (*offset) = i;
    
    return 1;
}

