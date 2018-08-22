/* pls check the header file: read_marker.h for more info about this function.                    */
/* Date 2013-04-26                                                                                */
#include        <fstream>
#include            <map>
#include         <vector>
#include        <cstring>
#include        <stdio.h>
#include       <stdlib.h>
#include        <ctype.h>
#include    "is_number.h"
#include "split_string.h"
#include      "globals.h"

bool read_marker(char* fmarker, unsigned long* num_given, unsigned long* num_inuse)
{
    if(fmarker==NULL)
    {
        printf("ERROR: Marker file is null, exited!\n");
        exit(1);
    }
    std::ifstream fileMarker (fmarker);
    if(!fileMarker.is_open())
    {
        printf("SNP file \'%s\' does NOT exist. Exited.\n", fmarker);
        exit(1);
    }
    if (verbose) 
    {
        ALLELE1.clear(); // ALT allele
        ALLELE2.clear(); // REF allele
        QUALITY1.clear();
        printf("Reading SNPs from file:\t\t\t%s...", fmarker);
    }
    
    char token='\t';                                                                     // caution!
    std::string line;
    while(fileMarker.good())
    {
        line.clear();
        getline(fileMarker, line);
        if(line.size() == 0)                        continue;
        if(line.find("#") != std::string::npos)   {printf("# annotation line skipped.\n");continue;}
        if(line.find("\t-\t") != std::string::npos) continue;                            // caution!

        (*num_given) ++;
        
        vector<std::string> infoline=split_string(line, token);
        if(infoline.size() < 5)
        {
            printf("SNP info line requires at least 5 columns: %s. Skipped!\n", line.c_str());
            continue;
        }

        std::string chr = infoline[1];
        if(CHR2SIZE.find(chr) == CHR2SIZE.end()) {continue;}
        std::string pos = infoline[2];
        if(!is_number(pos.c_str()))            
        {
            printf("Warning: position info of a SNP should be a number (read_marker(...)).\n"); 
            continue;
        }
        if (atol(pos.c_str()) > CHR2SIZE[chr]) 
        {   
            printf("Warning: position %ld of a SNP exceeds given chromosome size %ld. Skipped! \n", 
                   atol(pos.c_str()), CHR2SIZE[chr]);
            continue;
        }
        std::string pro = infoline[0]; // proj name
        std::string al1 = infoline[3]; // ref if background2 is not specified; or background allele
        std::string al2 = infoline[4]; // mut ..
        std::string qua = pro;         // 'proj name'
        if(infoline.size() >= 6) 
        {
            /* !!! caution: marker_score is a global variable                                     */
            if(atol(infoline[5].c_str()) < marker_score) continue;   // filter a marker with quality
            int addi = 5;
            while(addi < infoline.size())
            {
                qua  += ",";
                qua  += infoline[addi];                     // 'proj name,quality, cov, avg_hits...'
                addi ++;
            }
        }
        
        // note that alleles created by SHOREmap create (for outcross function) have been swapped in its output, if they are from parentB.
        // that is, such markers should be treated the same as those from parentA.
        
        /* record info about al1, al2 */
        if(background2 == 0)
        {
            ALLELE1.insert(std::pair<std::string, std::string>(chr+".#."+pos, al2)); // af to calculate
            ALLELE2.insert(std::pair<std::string, std::string>(chr+".#."+pos, al1));
        }
        else
        {
            ALLELE1.insert(std::pair<std::string, std::string>(chr+".#."+pos, al1)); // af to calculate
            ALLELE2.insert(std::pair<std::string, std::string>(chr+".#."+pos, al2));
        }
        /* QUALITY1 may only contain a project name if there is no quality info provided.         */
        // printf("%s\n", qua.c_str());
        QUALITY1.insert(std::pair<std::string, std::string>(chr+".#."+pos,qua));
        (*num_inuse) ++;
    }
    fileMarker.close();
    
    if(verbose) 
    {
        printf("done.\n");
        printf("Markers given to inuse:\t\t\t%ld: %ld.\n", *num_given, *num_inuse);
    }
    if (ALLELE1.size()>0 || ALLELE2.size()>0) 
         return  true; // at least one SNP
    else return false; // no SNP
}
