/* this function determines markers for the pool with SNPs of two parents against a reference seq. 
   Suppose parent a (pa) is with the same phenotype with F2, while parent b (pb) is not:               

   step 1: read in SNPs from parent a - collected in a map paSNP - done by func bg-ref-check;
   step 2: read in SNPs from parent b - collected in a map pbSNP - done by func bg-ref-check;
   step 3: remove common positions of paSNP and pbSNP;
           *note: can be improved by removing paSNP during reading and recording pbSNP

Note that not all of the SNPs in {parent_a + parent_b - parent_a*parent_b} will appear 
in F2 (that is, some parental SNP information is missed - number of markers in F2 is decreased). 
However, the following process guarantees all SNPs in parent_a and parent_b will be used as markers.

TODO: revise with extract consensus info of the pool for markers determined according to parents.

*/

#include        <assert.h>
#include          <time.h>
#include             <map>
#include          <vector>
#include         <fstream>
#include        <string.h>
#include        <iostream>
#include         <sstream>
#include        <stdlib.h>
#include       <algorithm>

#include        "globals.h"
#include      "is_number.h"
#include   "split_string.h"
#include  "check_ref_err.h"
#include    "read_marker.h"

using namespace std;

bool find_marker_pos_2parents_v3(char* paSNPfile, 
                                 char* pbREFfile, 
                                 char* pbSNPfile,
                                 char* paREFfile,
                                 char* F2SNPfile,
                                 unsigned long min_cov,
                                 unsigned long max_cov,
                                 double        min_ccd, 
                                 unsigned long max_N, 
                                 unsigned long max_ID, 
                                 unsigned long min_refq,
                                 std::string*  F2marker_ParentAsumBminusAB)
{
     /* 1.check SNPs of pa (with phenotype) with reference base call of pb                        */
     char*  fg_consencall_file = (char*)"";     // not used filtering step, set it as null file name
     std::string              pa_fg_snp_file_filtered = "";
     map<std::string, MARK6>  pa_fg_snp_map_filtered;
     pa_fg_snp_map_filtered.clear();
     if(paSNPfile != NULL)
     {    
         bgREFB4A.clear();
         check_ref_err(paSNPfile, 
                       fg_consencall_file, /* it is consensus-file related to paSNPfile - not used yet*/
                       pbREFfile,
                       min_cov,
                       max_cov,
                       min_ccd,
                       max_N, 
                       max_ID, 
                       min_refq, 
                       &pa_fg_snp_file_filtered, 
                       &pa_fg_snp_map_filtered);
         /* set flag as parent a                                                                      */
         map<std::string, MARK6>::iterator pa_snp_itr;
         map<std::string, MARK6>::iterator pa_snp_itr_end;
         pa_snp_itr     = pa_fg_snp_map_filtered.begin();
         pa_snp_itr_end = pa_fg_snp_map_filtered.end();
         while(pa_snp_itr != pa_snp_itr_end)
         {
             (*pa_snp_itr).second.par = "a";
             pa_snp_itr ++;
         }                       
         bgREFB4A.insert(bgREF.begin(), bgREF.end());
         bgREF.clear();         
         if(verbose)
         {
             cout << pa_fg_snp_file_filtered << " = bg-a markers supported by " << bgREFB4A.size();
             cout << " bg-ref-b." << endl;
         }
     }
    
    /* 2.check SNPs of pb (without phenotype) with reference base call of pa                      */
    /* TODO:  this can done combined with checking paSNP - now lazy to modify...                  */
    std::string              pb_fg_snp_file_filtered = "";
    map<std::string, MARK6>  pb_fg_snp_map_filtered;   
    pb_fg_snp_map_filtered.clear();
    if(pbSNPfile != NULL)
    {
        bgREFA4B.clear();
        check_ref_err(pbSNPfile,
                      fg_consencall_file, /* it is consensus-file related to pbSNPfile - not used yet*/
                      paREFfile, 
                      min_cov+1,
                      max_cov,
                      min_ccd,
                      max_N,
                      max_ID, 
                      min_refq,
                      &pb_fg_snp_file_filtered,
                      &pb_fg_snp_map_filtered);
        /* set flag as parent b                                                                   */
        map<std::string, MARK6>::iterator pb_snp_itr;
        map<std::string, MARK6>::iterator pb_snp_itr_end;
        pb_snp_itr     = pb_fg_snp_map_filtered.begin();
        pb_snp_itr_end = pb_fg_snp_map_filtered.end();
        while(pb_snp_itr != pb_snp_itr_end)
        {
            (*pb_snp_itr).second.par = "b";
            pb_snp_itr ++;
        }                  
        bgREFA4B.insert(bgREF.begin(), bgREF.end());
        bgREF.clear();                    
        if(verbose)
        {
            cout << pb_fg_snp_file_filtered << " = bg-b markers supported by " << bgREFA4B.size();
            cout << " bg-ref-a." << endl;
        }
    }
    cout << "Number of pa-markers - total: " << pa_fg_snp_map_filtered.size() << endl;
    cout << "Number of pb-markers - total: " << pb_fg_snp_map_filtered.size() << endl; 
    /* 3.remove common positions of paSNP and pbSNP                                              */
    if(pb_fg_snp_map_filtered.size() > 0 && pa_fg_snp_map_filtered.size() > 0)
    {
         map<std::string, MARK6>::iterator pbchp_itr_end;
         map<std::string, MARK6>::iterator   chp_itr;      // key=chr.#.pos of pa
         map<std::string, MARK6>::iterator   chp_itr_end;
         pbchp_itr_end  = pb_fg_snp_map_filtered.end();
           chp_itr      = pa_fg_snp_map_filtered.begin();
           chp_itr_end  = pa_fg_snp_map_filtered.end();
         while(chp_itr != chp_itr_end)
         {
             map<std::string, MARK6>::iterator pbchp_itr; // key=chr.#.pos of pb
             pbchp_itr = pb_fg_snp_map_filtered.find((*chp_itr).first);
             if(pbchp_itr != pbchp_itr_end)               // common discovered
             {
                 pb_fg_snp_map_filtered.erase(pbchp_itr);
                 if(!keep_common)
                 pa_fg_snp_map_filtered.erase(chp_itr++);
             }
             else
             {
                 chp_itr ++;
             }
        }
    }
    /* 4.check in case all SNPs of parents have been removed                                      */
    if(pb_fg_snp_map_filtered.size() == 0 && pa_fg_snp_map_filtered.size() == 0)
    {
        cout << "Warning: no markers available from given sets of markers of two parents!" << endl;
    }
    if(!keep_common)
    cout << "Number of pa-markers - specific (common removed): " << pa_fg_snp_map_filtered.size() << endl;
    else
    cout << "Number of pa-markers - specific (common kept)   : " << pa_fg_snp_map_filtered.size() << endl;    
    cout << "Number of pb-markers - specific (common removed): " << pb_fg_snp_map_filtered.size() << endl;
    /* 5.read and check F2-SNPs with paSNP and pbSNP, and record valid SNPs of F2 in file below   */
    /* merge two sets of snps */
    map<std::string, MARK6>  merged_fg_snp_map_filtered; // merge as one (so that elements are sorted)
    
    (*F2marker_ParentAsumBminusAB) = "";
    if(pa_fg_snp_map_filtered.size() > 0 && pb_fg_snp_map_filtered.size() > 0)
    {
        (*F2marker_ParentAsumBminusAB) = out_folder + "SHOREmap_created_F2Pab_specific.txt\0";
        merged_fg_snp_map_filtered.insert(pa_fg_snp_map_filtered.begin(), pa_fg_snp_map_filtered.end());
        pa_fg_snp_map_filtered.clear();
        merged_fg_snp_map_filtered.insert(pb_fg_snp_map_filtered.begin(), pb_fg_snp_map_filtered.end());
        pb_fg_snp_map_filtered.clear();    
    }
    else if(pa_fg_snp_map_filtered.size() > 0)
    {
        (*F2marker_ParentAsumBminusAB) = out_folder + "SHOREmap_created_F2onPa.txt\0";
        merged_fg_snp_map_filtered.insert(pa_fg_snp_map_filtered.begin(), pa_fg_snp_map_filtered.end());
        pa_fg_snp_map_filtered.clear();        
    }
    else if(pb_fg_snp_map_filtered.size() > 0)
    {
        (*F2marker_ParentAsumBminusAB) = out_folder + "SHOREmap_created_F2onPb.txt\0";
        merged_fg_snp_map_filtered.insert(pb_fg_snp_map_filtered.begin(), pb_fg_snp_map_filtered.end());
        pb_fg_snp_map_filtered.clear();        
    }
    /* sort the merged map according to chr and pos: merged_fg_snp_map_filtered;                  */
    map<std::string, map<unsigned long, MARK6> > twoParentSNP;
    
    map<std::string, MARK6>::iterator snp_itr;
    map<std::string, MARK6>::iterator snp_itr_end;
    snp_itr     = merged_fg_snp_map_filtered.begin();
    snp_itr_end = merged_fg_snp_map_filtered.end();
    while(snp_itr != snp_itr_end)
    {
        std::string chrpos = (*snp_itr).first;
        std::vector<std::string> keyinfo = split_string(chrpos, '#'); // chr.#.pos: caution: chr itself can have '.'; cannot split with '.'
        unsigned long            pos_tmp = atol(keyinfo[1].substr(1).c_str());
        
        map<std::string, map<unsigned long, MARK6> >::iterator twop_itr_tmp;
        twop_itr_tmp = twoParentSNP.find(keyinfo[0].substr(0, keyinfo[0].size()-1));
        
        if(twop_itr_tmp == twoParentSNP.end())
        {
            map<unsigned long, MARK6> map_tmp;
            map_tmp.insert(std::pair<unsigned long, MARK6>(pos_tmp, (*snp_itr).second));
            twoParentSNP.insert(std::pair<std::string, map<unsigned long, MARK6> >(keyinfo[0].substr(0, keyinfo[0].size()-1), map_tmp));
        }
        else
        {
            (*twop_itr_tmp).second.insert(std::pair<unsigned long, MARK6>(pos_tmp, (*snp_itr).second));
        }
        
        snp_itr ++;
    }
    
    /* 6.open file for writing SNPs of F2, which must pass checking by pa/b                       */
    fstream fout;
    fout.open ((*F2marker_ParentAsumBminusAB).c_str(), ios::out);
    if(!fout.is_open())
    {
        printf("Cannot open file to write filtered markers (in find_marker_pos_2parents(...)). ");
        printf("Exited. \n");
        exit(1);
    }
    unsigned long num_snp_created = 0;
    /* 7.open file for reading SNPs of F2 to check with pa/b                                      */
    /* get quality of base call in F2 at markers                                                  */
    bool ifF2marker  = true;
    FILE* fpt_marker = fopen(F2SNPfile, "r"); 
    if(fpt_marker   == NULL)
    {
        ifF2marker   = false;
    }
    unsigned long numgiven = 0;
    unsigned long numinuse = 0;
    marker_score = 0; // caution: no filtering
    if(ifF2marker && !read_marker(F2SNPfile, &numgiven, &numinuse))
    {
        cout << "ERROR: exited because no snps recorded from " << F2SNPfile << endl;
        exit(1);
    }
    if(ifF2marker && verbose)
    {
        cout << "output markers defined from parental lines with F2 quality information\n";
    }
    else
    {
        cout << "output markers without F2 quality information\n";
    }
    cout << "maximum quality score is found as " << quality_max << endl;
    map<std::string, map<unsigned long, MARK6> >::iterator sorted_snp_itr;
    map<std::string, map<unsigned long, MARK6> >::iterator sorted_snp_itr_end;
    sorted_snp_itr     = twoParentSNP.begin();
    sorted_snp_itr_end = twoParentSNP.end();
    while(sorted_snp_itr != sorted_snp_itr_end)
    {
        map<unsigned long, MARK6>::iterator sorted_pos_itr;
        map<unsigned long, MARK6>::iterator sorted_pos_itr_end;
        sorted_pos_itr     = (*sorted_snp_itr).second.begin();
        sorted_pos_itr_end = (*sorted_snp_itr).second.end();
        
        while(sorted_pos_itr != sorted_pos_itr_end)
        {
            MARK6     valmkr = (*sorted_pos_itr).second;
            
            /* as pa is with the same phenotype as F2, we observe its original allele frequency   */
            /* if change flag "FLAG4parentA/B" below, have to change at line137@read_allele_count2*/
            /* line119@print_filtered_marker, lines@plot_chr_winboost, lines@ShoreMap_outcross    
               line87@read_marker, line171@filter_with_marker_parent                              */            
            if(valmkr.par == "a")
            {
                fout << "FLAG4parentA" << "\t"; // project name: all the info below from parental lines
            }
            else
            if(valmkr.par == "b")
            {
                fout << "FLAG4parentB" << "\t"; // project name: all the info below from parental lines
            }
            else
            {
                cout << "Error: marker not proper." << endl;
                sorted_pos_itr ++;                   // //
                continue;
            }
            
            fout << (*sorted_snp_itr).first << "\t"; // chr
            fout << (*sorted_pos_itr).first << "\t"; // pos
            
            if(valmkr.par == "a")                    // 2013-04-23 17:07: caution! covs not swapped!
            {
                fout << valmkr.ref << "\t";          // ref base
                fout << valmkr.mut << "\t";          // alt base
            }
            else
            if(valmkr.par == "b")
            {
                fout << valmkr.mut << "\t";          // ref base
                fout << valmkr.ref << "\t";          // alt base
            }            
            
            string keystr("");
            keystr += (*sorted_snp_itr).first;
            keystr += ".#.";
            stringstream ss;
            ss.str("");
            ss << (*sorted_pos_itr).first;
            keystr += ss.str();
            
            if(ifF2marker)
            {
                map<string, string>::iterator qua_itr;
                qua_itr = QUALITY1.find(keystr);
                if(qua_itr != QUALITY1.end())
                {
                    vector<string> quainfo = split_string((*qua_itr).second, ',');
                    int iqua = 1;
                    while(iqua < quainfo.size())
                    {
                        fout << quainfo[iqua] << "\t";
                        iqua ++;
                    }
                }
                else                             // no qua, cov, af, avg_hit in quality_variant file
                {
                    int random_qua = (int)quality_max*0.6;
                    //srand(time(NULL));
                    int rnum = rand()%((int)quality_max-random_qua+1);
                    random_qua += rnum;
                    fout << random_qua << "\t20\t0\t1" << "\t";
                }
            }
            
            /* bg parent */
            fout << valmkr.qua << "\t";          // quality
            fout << valmkr.cov << "\t";          // coverage
            fout << valmkr.ccd << "\t";          // concordance
            if(!ifF2marker)
            fout << valmkr.aht << "\t";          // average hits
            
            if(ifF2marker)
            {
                map<std::string, std::string>::iterator bginfo_itr; 
                if(valmkr.par == "a")
                {
                    bginfo_itr = bgREFB4A.find(keystr);
                    assert(bginfo_itr != bgREFB4A.end());
                }
                else
                if(valmkr.par == "b")
                {
                    bginfo_itr = bgREFA4B.find(keystr);
                    assert(bginfo_itr != bgREFA4B.end());
                }
                string tmp("");
                tmp = (*bginfo_itr).second;
                std::replace(tmp.begin(), tmp.end(), '#', '\t');    
                fout << tmp << endl;
            }
            else
            {
                fout <<  endl;
            }

            num_snp_created ++;
            sorted_pos_itr ++;
        }        
        sorted_snp_itr ++;
    }
    fout.close();
    if (verbose) cout << num_snp_created << " markers have been created for the pool with 2 bgs." << endl;
    if(num_snp_created == 0) return false;
    else                     return true;
}

/*output format:
flag	        chr	pos	ref	mut	f2qua	f2cov	f2af	f2avg_hit  xqua	xcov	xaf	yqua	ycov	yaf

NOTE: the following information could be used to cluster F2 markers:
f2qua	(f2cov_from_consen_file)	xqua	xcov	xaf	yqua	ycov	yaf
*/
