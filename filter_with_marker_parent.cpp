/* pls check the header file: filter_with_marker_parent.h for more info about this function.      */
/* Date 2013-05-24                                                                                */
#include        <fstream>
#include            <map>
#include         <string>
#include         <vector>
#include        <stdio.h>
#include       <stdlib.h>

#include    "is_number.h"
#include "split_string.h"
#include      "globals.h"

using namespace std;

bool filter_with_marker_parent(char* fmarker_p, bool iphenotype)
{
    /* parent markers must satify the parameters to be used to keep/filter interested (F2)-markers:
       pmarker_score       =     0;
       pmarker_min_cov     =     0;
       pmarker_max_cov     =   INF;
       pmarker_min_freq    =   0.0;
    */
    unsigned long ipmarker_score; 
    unsigned long ipmarker_min_cov;
    unsigned long ipmarker_max_cov;
    double        ipmarker_min_freq;
    /* scale the cutoffs if necessary                                                             */
    if(!iphenotype && CMD.find("--pmarker-ab-ratio")!=CMD.end())
    {
        /* parent b using the scaled cutoffs from cmd line                                        */        
        double iscale[4]; // scaling base quality
        vector<std::string> infoscale = split_string(pmarker_ab_ratio, ',');
        if(infoscale.size() < 4)
        {
            printf("String %s for scaling of parent a/b not with size 4; last few ones set as 1.0\n",
                   (char*)pmarker_ab_ratio.c_str());
        }
        int iscl   = 0;
        while(iscl < 4)
        {
            bool r_from_cmd_valid = false; 
            if(iscl<infoscale.size())
            {
                if(is_number((char*)infoscale[iscl].c_str()))
                {
                    iscale[iscl] = atof(infoscale[iscl].c_str());
                    r_from_cmd_valid = true;
                }
                else
                {
                    printf("%s contains non-number characters, not used for scaling parent-ab. \n", 
                           (char*)infoscale[iscl].c_str());
                }
            }
            if(!r_from_cmd_valid)
            {
                iscale[iscl] = 1.0;
            }
	    iscl ++;
	}
        ipmarker_score   = (unsigned long)((double)pmarker_score*iscale[0]);
        ipmarker_min_cov = (unsigned long)((double)pmarker_min_cov*iscale[1]);
        ipmarker_max_cov = (unsigned long)((double)pmarker_max_cov*iscale[2]);
        ipmarker_min_freq= pmarker_min_freq*iscale[3];
    }
    else
    {
        /* parent a using the default cutoffs from cmd line                                       */
        ipmarker_score   = pmarker_score;
        ipmarker_min_cov = pmarker_min_cov;
        ipmarker_max_cov = pmarker_max_cov;
        ipmarker_min_freq= pmarker_min_freq;
    }

    unsigned long num_given = 0; // number of markers given in a file
    unsigned long num_inuse = 0; // number of markers to be used: kept or removed
    map<std::string, std::string> tmp_ALLELE1;
    map<std::string, std::string> tmp_ALLELE2;
    
    std::ifstream fileMarker_p (fmarker_p);
    if(!fileMarker_p.is_open())
    {
        printf("Marker file \'%s\' does NOT exist. Exited.\n", fmarker_p);
        exit(1);
    }
    if(verbose) printf("Filtering with parent-marker from file:\t\t%s...", fmarker_p);
    
    tmp_ALLELE1.clear();
    tmp_ALLELE2.clear();
    bool firstline = true;
    const char* token;
    std::string line;
    std::string mut_base;
    while(fileMarker_p.good())
    {
        line.clear();
        getline(fileMarker_p, line);
        if (line.length() == 0) {continue;}
        num_given ++;
        /* 0.projn 1.chr 2.pos 3.ref 4.mut 5.quality 6.raw_coverage 7.concordance 8.avg_hits      */
        vector<std::string> infoline = split_string(line, '\t');
        
        /* check necessary info-columns     */
        if(infoline.size() < 5)
        {
            printf("WARNING: at least 5 columns should be provided @line (if marker info)=%s.\n", 
                   (char*)line.c_str());
            continue;
        }
        /* if check min parent-marker score */
        if(ipmarker_score>0 && infoline.size()<6) 
        {
            printf("WARNING: at least 6 columns should be provided @line=%s.\n", 
                   (char*)line.c_str());
            //continue;
        }
        else
        if(atol(infoline[5].c_str()) < ipmarker_score)
        {
            /* set --marker-score as 0 (by default), if you dont want this. */
            continue;
        }
        /* if check min parent-marker coverage (raw) */
        if(CMD.find("--pmarker-min-cov") != CMD.end() || CMD.find("--pmarker-max-cov") != CMD.end())
        if(infoline.size() < 7)
        {
            printf("WARNING: at least 7 columns should be provided @line=%s.\n", 
            (char*)line.c_str());
        }
        else
        if(atol(infoline[6].c_str())<ipmarker_min_cov || atol(infoline[6].c_str())>ipmarker_max_cov)
        {
            /* set --pmarker-min-cov as 0   (by default), if you dont want this.
                                    & --pmarker-max-cov as INF (by default),
            */
            continue;
        }
        /* if check min parent-marker frequency/concordance */
        if(CMD.find("--pmarker-min-freq") != CMD.end())
        if(infoline.size() < 8)
        {
            printf("WARNING: at least 8 columns should be provided @line=%s.\n", 
            (char*)line.c_str());
        }
        else
        if(atof(infoline[7].c_str()) < ipmarker_min_freq)
        {
            /* set --pmarker-min-freq as 0.0 (by default), if you dont want this. */
            continue;
        }
        
        /* check chr id               */
        if(CHR2SIZE.find(infoline[1]) == CHR2SIZE.end()) {continue;}
        /* check position             */
        if (atol(infoline[2].c_str()) >  CHR2SIZE[infoline[1]]) {continue;}
        /* check if key exists        */
        std::string ale_id = infoline[1]+".#."+infoline[2];
        
        /* TODO: a parent-marker can be used if it safisfies quality, coverage, concordance       */
        
        /* if chrid and position are the same as the mutant alleles2 */
        map<std::string, std::string>::iterator ale_itr, ref_itr;
        ale_itr = ALLELE2.find(ale_id);
        if(ale_itr != ALLELE2.end())
        {
            ref_itr = ALLELE1.find(ale_id);
            /* if mut-base is the same    */
            if(background2) mut_base = infoline[3];
            else mut_base = infoline[4];
            if(mut_base == (*ale_itr).second)
            {
                if(iphenotype)
                {
                    /* keep   this marker common with parent a */
                    tmp_ALLELE1.insert(std::pair<string, string>((*ref_itr).first, (*ref_itr).second));
                    tmp_ALLELE2.insert(std::pair<string, string>((*ale_itr).first, (*ale_itr).second));
                }
                else
                {
                    /* remove this marker common with parent b */
                    ALLELE1.erase(ref_itr);
                    ALLELE2.erase(ale_itr);    
                }
                num_inuse ++; 
            }
        }
    }
    fileMarker_p.close();
    if(verbose) printf("done.\n");
    if(iphenotype)
    {
         /* keep   this marker */
         ALLELE1.clear();
         ALLELE2.clear();
         ALLELE1.insert(tmp_ALLELE1.begin(), tmp_ALLELE1.end());
         ALLELE2.insert(tmp_ALLELE2.begin(), tmp_ALLELE2.end());
         if(verbose) 
         {
              printf("Ratio of markers given to common-with-the-mutant (kept):\t%ld: %ld.\n", 
                     num_given, num_inuse);
         }
    }
    else
    {
         if(verbose) 
         {
              printf("Ratio of markers given to common-with-the-mutant (removed):\t%ld: %ld.\n", 
                     num_given, num_inuse);
         }
    }
    if (ALLELE1.size()>0 || ALLELE2.size()>0) 
         return  true; // at least one marker
    else return false; // no marker
}
