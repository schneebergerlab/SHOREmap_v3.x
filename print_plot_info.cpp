/* Function: print out the information that has been used in function: plot_chr_winboost for future
   independent manipulation.
   Date: 2012-12-23?
*/

#include  <stddef.h>
#include  <stdlib.h>
#include   <stdio.h>
#include   <fstream>
#include   <sstream>
#include   <cstring>
#include   <iomanip>

#include   "globals.h"

using namespace std;

int print_plot_info(std::string chrID,
                  double* sglPosiSet, double* sglFreqSet, unsigned long sglsize,
                  double* winPosiSet, double* winAFr2Set, double* winSTDnew2, double* winABSTSet, 
                  unsigned long awsize)
{
    bool appMode;
    /* single-marker frequency             */
    appMode = false;
    std::string sglstat = out_folder+"SHOREmap_stat_single_marker.txt\0";
    fstream fpout;
    std::ifstream ifile;
    ifile.open(sglstat.c_str());
    if (ifile.good())
    {
        printf("\tappending SHOREmap_stat_single_marker.txt with ");
        ifile.close();
        fpout.open(sglstat.c_str(), ios::out | ios::app);
        appMode = true;
    }
    else
    {
        printf("\tcreating SHOREmap_stat_single_marker.txt with ");
        fpout.open(sglstat.c_str(), ios::out);
    }
    if(!fpout.is_open())
    {
        printf("Cannot open file to write single-marker-stats (in print_plot_info.cpp). Exited.\n");
        exit(1);
    }
    else
    {
        printf("%ld individual markers\n", sglsize);    
        if(!appMode) fpout << "#chr\tpostion\tallele_frequency\n";  
        unsigned long num_sgl_entry = 0;
        while(num_sgl_entry < sglsize)
        {
            fpout << fixed << setprecision(0);
            fpout << chrID << "\t" << (*sglPosiSet) << "\t" << setprecision(8) << (*sglFreqSet) << endl;
            /* next entry */
            num_sgl_entry ++;
            sglPosiSet ++;
            sglFreqSet ++;
        }
        fpout.close();
    }
    
    /* window-marker frequency/boost-value */
    if(plot_window)
    {
        appMode = false;
        std::string winstat = out_folder+"SHROEmap_stat_window_markers.txt\0";
        fstream fpout2;
        ifile.open(winstat.c_str());
        if (ifile.good()) 
        {
            printf("\tappending SHROEmap_stat_window_markers.txt with ");
            ifile.close();
            fpout2.open(winstat.c_str(), ios::out | ios::app);
            appMode = true;
        }
        else
        {
            printf("\tcreating SHROEmap_stat_window_markers.txt with ");
            fpout2.open(winstat.c_str(), ios::out);
        }
        if(!fpout2.is_open())
        {
            printf("Cannot open file to write window-stats (in print_plot_info.cpp). Exited.\n");
            exit(1);
        }
        {
            printf("%ld windows of markers\n", awsize);
            if(!appMode) fpout2 << "#chr\tpostion\twin_allele_frequency\tCv\twin_boost_value\n";  
            unsigned long num_entry = 0;
            while(num_entry < awsize)
            {
                fpout2    << fixed << setprecision(0);
                fpout2    << chrID << "\t"   << *(winPosiSet+num_entry) << "\t";
                fpout2    << setprecision(8) << *(winAFr2Set+num_entry) << "\t";
                fpout2    << setprecision(8) << *(winSTDnew2+num_entry);
                if(plot_boost)
                {
                   fpout2 << "\t" <<  setprecision(8) << *(winABSTSet+num_entry) << endl;
                   ////winABSTSet ++;
                }
                else
                {
                   fpout2 << endl;
                }
                num_entry ++;
                /* next entry */
                ////winPosiSet ++;
                ////winAFr2Set ++;
            }
            fpout2.close();
        }
    }
    
    return 1;
}
