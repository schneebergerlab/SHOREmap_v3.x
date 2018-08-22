/* pls check the header file: print_filtered_marker.h for more info about this function.          */
/* Date 2013-05-28                                                                                */
#include     <fstream>
#include         <map>
#include      <vector>
#include     <cstring>
#include     <sstream>
#include     <stdio.h>
#include    <stdlib.h>

#include    "is_number.h"
#include "split_string.h"
#include      "globals.h"


bool print_filtered_marker(map<std::string, map<unsigned long, TRIPLE> > consen_info,
                           map<std::string, std::string> ref_allele,
                           map<std::string, std::string> mut_allele,
                           map<std::string, std::string> mut_quality)
{
    std::string marker_filter = "SHOREmap_marker_filtered_";
    /* note:   W: window; I: min; A: max; C: coverage; p: parent (a/b); F: allele frequency; 
             chr: chromosome; B: begin; E: end; N: number of markers of a window; S: base quality.
       E.G., WIN: minimum of number of markers required for a window(-analysis);
           chrIF: within given chr, the minimum allele frequency of a marker, etc.
    */
    if(CMD.find("--min-marker")               != CMD.end())     
        marker_filter += (".WIN"  + CMD["--min-marker"]);
    if(CMD.find("--min-coverage")             != CMD.end())     
        marker_filter += (".IC"   + CMD["--min-coverage"]);
    if(CMD.find("--max-coverage")             != CMD.end())     
        marker_filter += (".AC"   + CMD["--max-coverage"]);
    if(CMD.find("--marker-pa")                != CMD.end())     
        marker_filter += (".pa");
    if(CMD.find("--marker-pb")                != CMD.end())     
        marker_filter += (".pb");
    if(CMD.find("--pmarker-score")            != CMD.end())     
        marker_filter += (".pS"   + CMD["--pmarker-score"]);
    if(CMD.find("--pmarker-min-cov")          != CMD.end())  
        marker_filter += (".pIC"  + CMD["--pmarker-min-cov"]);
    if(CMD.find("--pmarker-max-cov")          != CMD.end())  
        marker_filter += (".pAC"  + CMD["--pmarker-max-cov"]);
    if(CMD.find("--pmarker-min-freq")         != CMD.end())  
        marker_filter += (".pIF"  + CMD["--pmarker-min-freq"]);
    if(CMD.find("--chromosome")               != CMD.end())        
        marker_filter += (".chr"  + CMD["--chromosome"]);
    if(CMD.find("--begin")                    != CMD.end())        
        marker_filter += (".chrB" + CMD["--begin"]);
    if(CMD.find("--end")                      != CMD.end())        
        marker_filter += (".chrE" + CMD["--end"]);
    if(CMD.find("--minfreq")                  != CMD.end())        
        marker_filter += (".chrIF"+ CMD["--minfreq"]);
    if(CMD.find("--maxfreq")                  != CMD.end())        
        marker_filter += (".chrAF"+ CMD["--maxfreq"]);
    if(CMD.find("--pmarker-ab-ratio")         != CMD.end())
        marker_filter += (".pabr" + CMD["--pmarker-ab-ratio"]);
    if(CMD.find("-bg-ref-filter")             != CMD.end())
        marker_filter += ".on_bgREF";
    else
        marker_filter += ".off_bgREF";
    marker_filter     += ".txt";
    marker_filter      = out_folder + marker_filter;
    printf("\tfiltered marker file:\t %s\n", marker_filter.c_str());
    fstream fout;
    fout.open (marker_filter.c_str(), ios::out);
    if(!fout.is_open())
    {
        printf("Cannot open file to write filtered markers (in print_filtered_marker.cpp). ");
        printf("Exited. \n");
        exit(1);
    }
    
    map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr;
    map<std::string, map<unsigned long, TRIPLE> >::iterator chr_itr_end;
    chr_itr     = consen_info.begin();
    chr_itr_end = consen_info.end();
    
    bool firstline = true;
    while(chr_itr != chr_itr_end)
    {
        map<unsigned long, TRIPLE>::iterator pos_itr;
        map<unsigned long, TRIPLE>::iterator pos_itr_end;
        pos_itr     = (*chr_itr).second.begin();
        pos_itr_end = (*chr_itr).second.end();
        while(pos_itr != pos_itr_end)
        {
            std::string key; // chr+".#."+pos
            std::stringstream ss;
            ss << (*pos_itr).first; // position: LONG to STRING type
            key = (*chr_itr).first + ".#." + ss.str();
            map<std::string, std::string>::iterator mkr_itr;
            mkr_itr = QUALITY1.find(key);
            if(mkr_itr == QUALITY1.end()) 
            {
                pos_itr ++;
                continue;
            }
            vector<std::string> projn_quality = split_string((*mkr_itr).second, ','); // projn_quality[0]=project name, projn_quality[1]=base_quality (if exists)
            
            if(firstline) firstline = false;
            else 
            fout << endl;
            fout << projn_quality[0]  << "\t"; // project name
            fout << (*chr_itr).first  << "\t"; // chromosome id
            fout << ss.str()          << "\t"; // position
            
            mkr_itr = ALLELE1.find(key);
            fout << (*mkr_itr).second << "\t"; // ref-base
            mkr_itr = ALLELE2.find(key);
            fout << (*mkr_itr).second << "\t"; // mut-base
            if(projn_quality.size()>1)
            fout << projn_quality[1]  << "\t"; // quality
            
            //(*pos_itr).second.Ci[0] ; ref_count
            //(*pos_itr).second.Ci[1] ; mut_count
            //(*pos_itr).second.Ci[2] ; err_count
            unsigned long cov_refplusmut = (*pos_itr).second.Ci[0]+(*pos_itr).second.Ci[1];
            
            if(QUALITY1[key].find("FLAG4parentB") != std::string::npos)
            {
                /* 2013-09-20: from background parent b, which does not have the phenotype        */
                /* in this case, when reading consensus info, we swapped the count the ref and mut*/
                /* to show allele freq of one bg; for record and check, we swap it back           */
                fout << (*pos_itr).second.Ci[0] << "\t"; // coverage on mut -- caution of only mut
                fout << (double)(*pos_itr).second.Ci[0]/cov_refplusmut << "\t";
                
                // 2013-07-09 - no above-swap when recording consensus info, thus normal output
                // fout << (*pos_itr).second.Ci[1] << "\t"; // coverage on mut
                // fout << (double)(*pos_itr).second.Ci[1]/cov_refplusmut << "\t";
            }
            else
            {  
                /* normal case: from background parent a, which has the same phenotype as mutant  */        
                fout << (*pos_itr).second.Ci[1] << "\t"; // coverage on mut -- caution of only mut
                fout << (double)(*pos_itr).second.Ci[1]/cov_refplusmut << "\t";
            }
            fout << "1";                                 // avg_hit 
            pos_itr ++;
        }
        chr_itr ++;
    }

    fout.close();
    return true;
}
