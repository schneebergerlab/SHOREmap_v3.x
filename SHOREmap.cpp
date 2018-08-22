////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  SHOREMap version 3.x: for fast and accurate identification of causal mutations in plants      //
//  Written by:           Hequan Sun, Korbinian Schneeberger                                      //
//  Based on:             version 2.0 in perl and R, and version 2.1,3.x written in C/C++         //
//  Date:                 2013-03-01 to 2017-09-26                                                //
//                                                                                                //
//  New features:         TO-SPECIFY                                                              //
//  SHOREmap is available under GPL license                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////COPYRIGHT//////////////////////////////////////////////////
/* SHOREMap version 3.x: for fast and accurate identification of causal mutations in plants 
    Copyright (C) 2013~2017  <Hequan Sun, Korbinian Schneeberger>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/
*/
////////////////////////////////////////////////////////////////////////////////////////////////////

#include             <stddef.h>
#include             <stdlib.h>
#include              <stdio.h>
#include               <math.h>
#include               <string>
#include             <string.h>
#include              <ctype.h>
#include               <time.h>
#include            "globals.h"
#include      "init_outcross.h"
#include   "ShoreMap_extract.h"
#include  "ShoreMap_outcross.h"
#include     "init_backcross.h"
#include  "ShoreMap_annotate.h"
#include   "ShoreMap_create2.h" //2013-09-17
#include   "ShoreMap_convert.h"
#include  "ShoreMap_idFilter.h"

using namespace std;

// declare sub functions
void init_global(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    double startT  = clock();
    init_global(argc, argv);
    double finishT = clock();
    
    printf("\nTime consumed %.4f seconds.\n", (finishT-startT)/1000000.0);
    return 0;
}
// define functions
void init_global(int argc, char* argv[])
{
    // this function read command line parameters
    string version = "3.6";
    std::string usage_global = "";
    usage_global += "Usage of SHOREmap version " + version + ": SHOREmap FUNCTION [OPTIONS]\n\n";
    usage_global += " FUNCTION is one of:\n\n";
    //usage_global += "  outcross \tfor analysis of outcross-mapping population\n";
    usage_global += "  backcross\tanalyzes  backcrossing-derived mapping population\n";
    usage_global += "  outcross \tanalyzes  outcrossing--derived mapping population\n";
    usage_global += "  annotate \tannotates mutations based  on  gene    feature \n";
    usage_global += "  create   \tcreates   markers   according  to  two backgrounds \n";
    usage_global += "  extract  \textracts  base call info     according to markers \n";
    usage_global += "  convert  \tconverts  VCF4.1   into     acceptable formats; or \n";
    usage_global += "           \tconverts  shore pileup into shore consensus_summary format \n";
    usage_global += "  idFilter \tfilters   indels in one sample according to another\n";
    usage_global += "  about    \tdisplays  general/version/publication  info about SHOREmap\n";
    usage_global += "\n\nFor detailed OPTIONS, try: SHOREmap FUNCTION.\n\n";
    std::string about = "";
    about += "\nGeneral: development of SHOREmap was started in the lab of Detlef Weigel, and\n";
    about += "  is also continued in the lab of Korbinian Schneeberger.\n\n";
    about += "Latest version: " + version +", developed by Hequan Sun and Korbinian Schneeberger\n\n";
    about += "Publications on SHOREmap:\n";
    about += "  *Version 3.0 with applications is introduced here: \n";
    about += "    H.Q. Sun and K. Schneeberger (2015)\n";
    about += "    SHOREmap v3.0: fast and accurate identification of causal mutations from forward\n";
    about += "    genetic screens. In: Alonso, J.M., and Stepanova, A.N. (Eds.), \n";
    about += "    Plant Functional Genomics: Methods and Protocols (pp381-395), 2nd Edition, New York: Springer.\n";    
    about += "  *Original:\n";
    about += "    K. Schneeberger, S. Ossowski, C. Lanz, T. Juul, et al. (2009)\n";
    about += "    SHOREmap: simultaneous mapping and mutation identification by deep sequencing.\n";
    about += "    Nature Methods. 6, 550-551.\n";
    about += "  *Analysis of outcrossing data is first introduced here\n";
    about += "    V.C. Galvão, K.J. Nordström, C. Lanz, P. Sulz, et al. (2012) \n";
    about += "    Synteny-based Mapping-by-Sequencing enabled by Targeted Enrichment. Plant J.\n";
    about += "  *Analysis of backcrossing data is described here:\n";
    about += "    B. Hartwig, G.V. James, K. Konrad, K. Schneeberger, F. Turck et al. (2012) \n";
    about += "    Fast isogenic mapping-by-sequencing of EMS-induced mutant bulks. Plant Physiology.\n";

    if (argc < 2) 
    {
        printf("\nNOT enough input parameters. Check necessary parameters below: \n\n");
        printf("%s\n", usage_global.c_str());
        exit(1);
    }
    
    char* myTask = argv[1];
    if(!strcmp(myTask, "backcross"))
    {
        if(argc > 2)
        printf("\nCalling %s: \n",
         myTask);
        strcatCMD += "backcross ";
        init_backcross(argc, argv);
        printf("\nBackcross analysis done.\n");
    }
    else if(!strcmp(myTask, "outcross"))
    {
        if(argc > 2)
        printf("\nCalling %s: \n", myTask);
        strcatCMD += "outcross ";
        if(ShoreMap_outcross(argc, argv))
        printf("\nOutcross analysis done.\n");
    }
    else if(!strcmp(myTask, "annotate"))
    {
        if(argc > 2)
        printf("\nAnnotating markers (SNPs): \n");
        strcatCMD += "annotate ";
        ShoreMap_annotate(argc, argv);
        printf("\nAnnotation done.\n");
    }
    else if(!strcmp(myTask, "convert"))
    {
        if(argc > 2)
        printf("\nConverting format: \n\n");
        strcatCMD += "convert ";
        ShoreMap_convert(argc, argv);
        printf("\nConvertion done.\n");
    }
    else if(!strcmp(myTask, "create"))
    {
        if(argc > 2)
        printf("\nCreating F2-markers based on background lines: \n");
        strcatCMD += "create ";
        ShoreMap_create2(argc, argv);                          //2013-09-17
        printf("\nCreation of markers done.\n");
    }
    else if(!strcmp(myTask, "extract"))
    {
        if(argc > 2)
        printf("\nExtracting consensus according to markers: \n");
        strcatCMD += "extract ";
        ShoreMap_extract(argc, argv);
        printf("\nEXtraction done.\n");
    }
    else if(!strcmp(myTask, "idFilter"))
    {
        if(argc > 2)
        printf("\nFiltering indels in one file according to another: \n");
        strcatCMD += "idFilter ";
        ShoreMap_idFilter(argc, argv);
        printf("\nIndel-filering done.\n");
    }
    else if(!strcmp(myTask, "about"))
    {
        strcatCMD += "about\n";
        printf("%s\n", about.c_str());
        exit(0);
    }
    else
    {
        printf("\nINVALID command: %s. Check COMMANDs below: \n\n", myTask);
        printf("%s\n", usage_global.c_str());
    }
}
