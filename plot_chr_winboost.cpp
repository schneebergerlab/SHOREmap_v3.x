/* this function plots single-marker/window-marker (boost) allele frequency along a chromosome    */
/* Date  : 2013-Apr-24                                                                            */
/* Date  : 2014-Feb-20 : test    3D-plot            with     quality      score     on    markers */
/* Date  : 2014-Mar-18 : test    k-means            for                   ranking         markers */
/* Date  : 2014-May-05 : use af=alt/(alt+ref+error) for      winAF2new                            */
/* Date  : 2014-Sep-05 : only pull down boost-values of  a  chromosome  when  no  confident  peak */
/* Date  : 2015-Jun-15 : add     PCI     plotting   for      backcrossing                    data */
#include                      <stddef.h>
#include                      <stdlib.h>
#include                       <stdio.h>
#include                        <math.h>
#include                        <string>
#include                      <string.h>
#include                      <assert.h>
#include                           <map>
#include                        <vector>
#include                       <fstream>
#include                       <sstream>
#include                      <iostream>
#include                     "globals.h"
#include             "print_plot_info.h"
#include                "split_string.h"
#include                      "kmeans.h"
#include                "find_z_value.h"
#include  "Calc_Wilson_score_interval.h"
// visualization
#include             "dislin/dislin_d.h"

using namespace std;

struct WININFO
{
    double pos;                                                              // window        center
    double val;                                                              // window   average/std
    double beg;                                                              // window         begin
    double end;                                                              // window           end
};

bool plot_chr_winboost(std::string                    chrID, 
                       unsigned long           plot_ith_chr,
                       map<unsigned long, TRIPLE> interData,
                       map<unsigned long, TRIPLE>  filtered)
{
    double  interest_regPos[2];
    double  interest_regHET[2];
    char    iregion[1024];
    double  avg_threshold = interval_min_mean;
    double   cv_threshold = interval_max_cvar;
    //
    double  max_quality = 0.0;
    double  miax_bg1_cov[2];
    double  cluster_avg_cov_bg1;
    double  cluster_dim_cov_bg1;
    double  miax_bg2_cov[2];
    double  cluster_avg_cov_bg2;
    double  cluster_dim_cov_bg2;   
    //
    unsigned long  dtSize = interData.size();
    double* myPosiSet = (double*)malloc((dtSize+1)*sizeof(double));
    double* myFreqSet = (double*)malloc((dtSize+1)*sizeof(double));
    double* myScorSet = (double*)malloc((dtSize+1)*sizeof(double));
    
    // PCI plot - 2015-Jun-15 17:52
    double* PciPosiSet;
    double* PciFreqSet;
    double* PciEuppSet;
    double* PciElowSet;
    double  oberveredAFMean = 0.0;
    if(pci && pci_chr.compare(chrID)==0)
    {
        PciPosiSet   = (double*)malloc((dtSize+1)*sizeof(double));
        PciFreqSet   = (double*)malloc((dtSize+1)*sizeof(double));
        PciEuppSet   = (double*)malloc((dtSize+1)*sizeof(double));
        PciElowSet   = (double*)malloc((dtSize+1)*sizeof(double));
        percentile_z = find_z_value(z_table, pci_cfd);
    }
    
    /* prepare for clustering of markers..........................................................*/
    unsigned long       n = dtSize;
    unsigned long       m = 0;
    unsigned long       k = (clusterK>dtSize)?dtSize:clusterK;
    double              t = 0.01;
    double** markerData;    
    double**  centroids;
    unsigned long* cluster_counts;      
    assert(myPosiSet && myFreqSet && myScorSet);
    // test condition: when k>=2 ...................................................................
    if(k >= 2)
    if(strcatCMD.find("SHOREmap outcross") != std::string::npos)
    {
        std::vector<string> quainfo = split_string((*(QUALITY1.begin())).second, ',');
        if(quainfo.size() == 11)            // caution! generated format by SHOREmap create function
        {
            m = 8;
        }
        else if(quainfo.size() <= 5)        // caution! generated format by  SHORE  variant  calling
        {
            m = 2;
        }
    }
    else if(strcatCMD.find("SHOREmap backcross")  != std::string::npos)
    {
        if(!bg_ref_filter && strcatCMD.find("--bg")==std::string::npos)
        {
            /* clustering..with..foreground..quality..and.....coverage.....of.....mutant.....base */
            m = 2;
        }
        else if(!bg_ref_filter && strcatCMD.find("--bg")!= std::string::npos)
        {
            /* clutering..with..foreground+background.quality and..coverage...of...mutant....base */
            m = 4;
        }
        else if(bg_ref_filter)
        {
            /* clustering with fg quality, coverage of mutbase, bg quality and cov, af of refbase */
            m = 5;
        }
    }
    else 
    {
        printf("ERROR: dimension of markerData is unknown??? Exited.\n");
        exit(1);
    }
    if(k >= 2)
    {
        markerData     = (double**)malloc((dtSize+1)*sizeof(double*));
        centroids      = (double**)malloc((k+1)*sizeof(double*));
        cluster_counts = (unsigned long*)malloc((k+1)*sizeof(unsigned long));
        assert(markerData && centroids && cluster_counts);
        for(int mali = 0; mali < dtSize+1; mali ++)
        {
            markerData[mali]    = (double*)malloc(m*sizeof(double));
            assert(markerData[mali]);
            if(mali <= k)
            {
                centroids[mali] = (double*)malloc(m*sizeof(double));
                assert(centroids[mali]);
            }
        }
    }
    /* prepare...xray:position....and.....yray....(and zray):....frequency...for.........plotting */
    miax_bg1_cov[0] = INF;
    miax_bg1_cov[1] = 0;
    miax_bg2_cov[0] = INF;
    miax_bg2_cov[1] = 0;
    map<unsigned long, TRIPLE>::iterator mkr_itr     = interData.begin();
    map<unsigned long, TRIPLE>::iterator mkr_itr_end = interData.end();
    unsigned long ipos    = 0;
    unsigned long ipci    = 0;
    while(mkr_itr != mkr_itr_end)
    {
        /* this (position, frequency): caution with error ........................................*/
        double icov       = (double)(*mkr_itr).second.Ci[0]+(*mkr_itr).second.Ci[1]+(*mkr_itr).second.Ci[2];
        *(myPosiSet+ipos) = (double)(*mkr_itr).first;
        *(myFreqSet+ipos) = (double)(*mkr_itr).second.Ci[1]/icov;           // ref/icov or mut/icov 
        
        if(pci && pci_chr.compare(chrID)==0)
        {
            if((*mkr_itr).first>=pci_start && (*mkr_itr).first<=pci_end)
            {
                *(PciPosiSet+ipci) = (double)(*mkr_itr).first;
                *(PciFreqSet+ipci) = (double)(*mkr_itr).second.Ci[1]/icov;
                double interval[2];
                Calc_Wilson_score_interval(icov, (*mkr_itr).second.Ci[1], percentile_z, interval);
                *(PciEuppSet+ipci) = interval[1] - *(PciFreqSet+ipci); 
                *(PciElowSet+ipci) = *(PciFreqSet+ipci) - interval[0];
                oberveredAFMean   += *(PciFreqSet+ipci);
                ipci ++;
            }
        }
        
        // new.....2014-02-20: find chr+".#."+pos and assign quality/rank score ..................*/
        // m = 2;
        if(k >= 2)
        {
            markerData[ipos][0] = (double)(*mkr_itr).second.Ci[1];          // fg-mut-cov...........
            std::stringstream posTmp;
            posTmp << (unsigned long)(*mkr_itr).first;
            string qkey("");
            qkey += chrID;
            qkey += ".#.";
            qkey += posTmp.str();
            map<std::string, std::string>::iterator qitr = QUALITY1.find(qkey); 
            string qua = (*qitr).second.c_str();                            // like qua = eop002,40.
            std::vector<string> qinfo = split_string(qua, ',');
            markerData[ipos][1]       = atof(qinfo[1].c_str());             // fg-mut-qua...........
            if(markerData[ipos][1] > max_quality)
            {
                max_quality = markerData[ipos][1];
            }
           
            if(strcatCMD.find("SHOREmap outcross") != std::string::npos     && m==8)
            {        
                // qinfo <- FLAG4parentA,(40,20,0.571429,1),(40,36,1),(40,29,1).....................    
                markerData[ipos][2] = atof(qinfo[5].c_str());  // qua
                markerData[ipos][3] = atof(qinfo[6].c_str());  // cov
                markerData[ipos][4] = atof(qinfo[7].c_str());  // af
                markerData[ipos][5] = atof(qinfo[8].c_str());  // qua
                markerData[ipos][6] = atof(qinfo[9].c_str());  // cov
                markerData[ipos][7] = atof(qinfo[10].c_str()); // af
                
                if(miax_bg1_cov[0] > markerData[ipos][3]) miax_bg1_cov[0] = markerData[ipos][3];
                if(miax_bg1_cov[1] < markerData[ipos][3]) miax_bg1_cov[1] = markerData[ipos][3];
                
                if(miax_bg2_cov[0] > markerData[ipos][6]) miax_bg2_cov[0] = markerData[ipos][6];
                if(miax_bg2_cov[1] < markerData[ipos][6]) miax_bg2_cov[1] = markerData[ipos][6];
            }
            else if(strcatCMD.find("SHOREmap backcross")!=std::string::npos && m==4)
            {
                /*fg-cov, fg-qua, bg-cov, bg-cov................................................. */
                map<std::string, std::string>::iterator bgmkr_itr = bgMARKER.find(qkey);
                std::vector<string> bgmkr_info = split_string((*bgmkr_itr).second, '#');
                markerData[ipos][2] = atof(bgmkr_info[1].c_str());                     // bg-mut-cov
                markerData[ipos][3] = atof(bgmkr_info[0].c_str());                     // bg-mut-sco
            }
            else if (strcatCMD.find("SHOREmap backcross")!=std::string::npos && m==5)
            {
                /* clustering with fg quality, cov of mutbase, bg quality and cov, af of ref base */
                map<std::string, std::string>::iterator bgref_itr = bgREF.find(qkey);
                std::vector<string> bgref_info = split_string((*bgref_itr).second, '#');
                markerData[ipos][2] = atof(bgref_info[1].c_str());                     // bg-ref-cov
                markerData[ipos][3] = atof(bgref_info[0].c_str());                     // bg-ref-sco
                markerData[ipos][4] = atof(bgref_info[2].c_str());                     // bg-ref-ccd
            }
            else 
            {
                ;
            } 
        }       
        /* next..........................................................................position */
        ipos ++;
        mkr_itr ++;
    }
    oberveredAFMean = oberveredAFMean/ipci;
    /* scale the factors: column 0 and 2 (if bg-mut): (cov-mean)/mean, column 1 and 3: quality/40 */
    //fg-cov, fg-qua
    if(verbose && k>=2)
    {
        cout << "\t" << m << " items will be used for clustering markers with " << endl;
        cout << "\t\tfgmut_avg_cov=" << cluster_avg_coverage << ", ";
        cout <<     "fgmut_dim_cov=" << cluster_dim_coverage << endl;
    }
    double bgref_cdim;  // caution: further accessed when re-scaling
    double bgref_cmean;
    unsigned long oni;
    if(k >= 2)
    if(m == 2)
    {
        oni = 0;
        while(oni < interData.size())
        {
            ////printf("%ld: %.0f %.0f --> ",  oni, markerData[oni][0],  markerData[oni][1]);
            markerData[oni][0] = 1-0.5*fabs(markerData[oni][0]-(double)cluster_avg_coverage)/(double)cluster_dim_coverage;
            markerData[oni][1] = markerData[oni][1]/max_quality;
            ////printf("%.2f %.0f \n",  markerData[oni][0],  markerData[oni][1]);
            oni++;
        }
    }
    else if(strcatCMD.find("SHOREmap outcross")!=std::string::npos && m==8)
    {
        assert(miax_bg1_cov[0] <=  miax_bg1_cov[1] && miax_bg1_cov[0]>0);
        cluster_avg_cov_bg1     = (miax_bg1_cov[1] +  miax_bg1_cov[0])/2.0;
        cluster_dim_cov_bg1     = (miax_bg1_cov[1] -  miax_bg1_cov[0])/2.0;
        if(cluster_dim_cov_bg1 == 0) cluster_dim_cov_bg1 = 1;   // caution!
        assert(miax_bg2_cov[0] <=  miax_bg2_cov[1] && miax_bg2_cov[0]>0);
        cluster_avg_cov_bg2     = (miax_bg2_cov[1] +  miax_bg2_cov[0])/2.0;
        cluster_dim_cov_bg2     = (miax_bg2_cov[1] -  miax_bg2_cov[0])/2.0;
        if(cluster_dim_cov_bg2 == 0) cluster_dim_cov_bg2 = 1;   // caution!
        if(verbose)
        {
            cout << "\t\tbg1_avg_cov=" << cluster_avg_cov_bg1 << ", ";
            cout <<     "bg1_dim_cov=" << cluster_dim_cov_bg1 << endl;
            cout << "\t\tbg2_avg_cov=" << cluster_avg_cov_bg2 << ", ";
            cout <<     "bg2_dim_cov=" << cluster_dim_cov_bg2 << endl;
        }
        oni = 0;
        while(oni < interData.size())
        {
            ////printf("\t%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f -> ",
            ////       markerData[oni][0], markerData[oni][1],  
            ////       markerData[oni][2], markerData[oni][3], markerData[oni][4],  
            ////       markerData[oni][5], markerData[oni][6], markerData[oni][7]);
            
            markerData[oni][0] = 1-0.5*fabs(markerData[oni][0]-(double)cluster_avg_coverage)/(double)cluster_dim_coverage;
            markerData[oni][1] = markerData[oni][1]/max_quality;
            markerData[oni][2] = markerData[oni][2]/max_quality;
            markerData[oni][3] = 1-0.5*fabs(markerData[oni][3]-(double)cluster_avg_cov_bg1)/(double)cluster_dim_cov_bg1;
            // markerData[ipos][4] = markerData[ipos][4];
            markerData[oni][5] = markerData[oni][5]/max_quality;
            markerData[oni][6] = 1-0.5*fabs(markerData[oni][6]-(double)cluster_avg_cov_bg2)/(double)cluster_dim_cov_bg2;
            // markerData[ipos][7] = markerData[ipos][7];
            
            ////printf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n",  
            ////       markerData[oni][0], markerData[oni][1],  
            ////       markerData[oni][2], markerData[oni][3], markerData[oni][4],  
            ////       markerData[oni][5], markerData[oni][6], markerData[oni][7]);
            oni ++;
        }
    }
    else if(strcatCMD.find("SHOREmap backcross")!=std::string::npos && m>=4 && m<=5)
    {
        // fg-cov, fg-qua, bg-cov, bg-cov
        if(m == 5)
        {
            if(strcatCMD.find("--bg-ref-cov-max") == std::string::npos &&
               strcatCMD.find("--bg-ref-cov") == std::string::npos)
            {
                printf("Warning: no coverage constraints on bg-ref bases, default will apply!\n");
            }
            else if(strcatCMD.find("--bg-ref-cov-max") == std::string::npos &&
                    strcatCMD.find("--bg-ref-cov") != std::string::npos)
            {
                    printf("Warning: only constraint on min-cov of ref base, default will apply!\n");
            }
            else if (strcatCMD.find("--bg-ref-cov-max")!=std::string::npos &&
                     strcatCMD.find("--bg-ref-cov")    ==std::string::npos)
            {
                    printf("Warning: no min-cove constraint on bg-ref bases, default will apply!\n");
            }
            else ;
            bgref_cdim  = (bg_ref_cov_max - bg_ref_cov)/2.0;
            bgref_cmean = (bg_ref_cov_max + bg_ref_cov)/2.0;
            if(bgref_cdim == 0.0) bgref_cdim = 1.0;//caution!
            if(verbose)
            printf("\tbgref_avg_coverage=%.0f, bgref_dim_coverage=%.0f\n", bgref_cmean, bgref_cdim);
        }
        oni = 0;
        while(oni < interData.size())
        {
            ////printf("\t%.0f %.0f %.0f %.0f",
            ////       markerData[oni][0], markerData[oni][1],  
            ////       markerData[oni][2], markerData[oni][3]);
            ////if(m == 5)
            ////{
            ////    printf(" %.2f\t\t -->\t", markerData[oni][4]);
            ////}
            ////else
            ////{
            ////    printf("\t\t -->\t");
            ////}
            markerData[oni][0] = 1-0.5*fabs(markerData[oni][0]-(double)cluster_avg_coverage)/(double)cluster_dim_coverage;
            markerData[oni][1] = markerData[oni][1]/max_quality;   
            if(m == 4)     
            {
                markerData[oni][2] = 1-markerData[oni][2]/(double)bg_read;                          // lower  - better   
                markerData[oni][3] = (max_quality-markerData[oni][3])/max_quality;                  // bg-mut/lower  - better  
            }
            else // if (m == 5)
            {
                markerData[oni][2] = 1-0.5*fabs(markerData[oni][2]-bgref_cmean)/bgref_cdim;          // higher - better
                markerData[oni][3] = markerData[oni][3]/max_quality;                                 // bg-ref-higher  - better
            }
            ////printf("%.2f %.2f %.2f %.2f",  
            ////       markerData[oni][0], markerData[oni][1], 
            ////       markerData[oni][2], markerData[oni][3]);
            ////if(m == 5)
            ////{
            ////    printf(" %.2f\n", markerData[oni][4]);
            ////}
            ////else
            ////{
            ////    printf("\n");
            ////}
            oni ++;
        }
    }
    else
    {
        ;
    }   
    char** real_value_centeroids;
    if(k >= 2)
    {
        if(verbose) printf("\tScaling of factors for %ld markers done.\n", dtSize);
        real_value_centeroids = (char**)malloc((k+1)*sizeof(char*));
        assert(real_value_centeroids);
        for(int irvc=0; irvc<k+1; irvc++)
        {
            real_value_centeroids[irvc] = (char*)malloc(256*sizeof(char));//caution!
            assert(real_value_centeroids[irvc]);
        }
        k_means(markerData, n, m, k, t, centroids, cluster_counts, myScorSet);
        if(verbose) 
        printf("\t%ld markers ==> %ld (actual; %ld expected) clusters. ", dtSize, k, clusterK);
        
        
        /* recover re-scaled values of centroids to real ones.....................................*/
        if(strcatCMD.find("SHOREmap outcross")!=std::string::npos && (m==2 || m==8))
        {
            if(verbose) printf("Non-scaled centers of clusters given in PDF.\n");
            double cen_fg_cov_low, cen_fg_cov_high, cen_fg_qua;
            double cen_bg1_cov_low, cen_bg1_cov_high, cen_bg1_cov, cen_bg1_qua;
            double cen_bg2_cov_low, cen_bg2_cov_high, cen_bg2_cov, cen_bg2_qua;
            /*
            markerData[oni][0] => fg-cov
            markerData[oni][1] => fg-qua  
            markerData[oni][2] => bg1-qua 
            markerData[oni][3] => bg1-cov
            markerData[ipos][4]=> bg1-af
            markerData[oni][5] => bg2-qua
            markerData[oni][6] => bg2-cov
            markerData[ipos][7]=> bg2-af
            */
            for(int ki=k-1; ki>=0; ki --)
            {
                printf("\t%.2f %.2f ",
                       centroids[ki][0],
                       centroids[ki][1]);
                if(m == 8)
                printf("%.2f %.2f %.2f %.2f %.2f %.2f -> ",
                       centroids[ki][2],
                       centroids[ki][3],
                       centroids[ki][4],
                       centroids[ki][5],
                       centroids[ki][6],
                       centroids[ki][7]);
                else
                printf(" -> ");
                
                cen_fg_cov_low      = (1-centroids[ki][0])*(double)cluster_dim_coverage;
                cen_fg_cov_low      = cen_fg_cov_low*2.0;
                cen_fg_cov_high     = (double)cluster_avg_coverage + cen_fg_cov_low;
                cen_fg_cov_low      = (double)cluster_avg_coverage - cen_fg_cov_low;
                if(cen_fg_cov_low   < (double)filter_min_coverage)
                {
                    cen_fg_cov_low  = cen_fg_cov_high;
                }
                cen_fg_qua          = centroids[ki][1]*max_quality;
                
                if(m == 8)
                {
                    cen_bg1_qua         = centroids[ki][2]*max_quality;
                    cen_bg1_cov_low     = (1-centroids[ki][3])*(double)cluster_dim_cov_bg1;
                    cen_bg1_cov_low     = cen_bg1_cov_low*2.0;
                    cen_bg1_cov_high    = (double)cluster_avg_cov_bg1 + cen_bg1_cov_low;
                    cen_bg1_cov_low     = (double)cluster_avg_cov_bg1 - cen_bg1_cov_low;
                    if(cen_bg1_cov_low  < (double)miax_bg1_cov[0])
                    {
                        cen_bg1_cov_low = cen_bg1_cov_high;
                    }
                
                    cen_bg2_qua         = centroids[ki][5]*max_quality;
                    cen_bg2_cov_low     = (1-centroids[ki][6])*(double)cluster_dim_cov_bg2;
                    cen_bg2_cov_low     = cen_bg2_cov_low*2.0;
                    cen_bg2_cov_high    = (double)cluster_avg_cov_bg2 + cen_bg2_cov_low;
                    cen_bg2_cov_low     = (double)cluster_avg_cov_bg2 - cen_bg2_cov_low;
                    if(cen_bg2_cov_low  < (double)miax_bg2_cov[0])
                    {
                        cen_bg2_cov_low = cen_bg2_cov_high;
                    }
                }
                ////
                if(cluster_counts[ki] > 0)
                {
                    sprintf(real_value_centeroids[ki], 
                            "    cluster %2d: (%-3.1f, %3.1f) %2.1f", ki,
                            cen_fg_cov_low, cen_fg_cov_high, cen_fg_qua);
                    if(m == 8)
                    sprintf(real_value_centeroids[ki], 
                            "%s %2.1f (%-3.1f, %3.1f) %.2f %2.1f (%-3.1f, %3.1f) %.2f %8ld\0", 
                            real_value_centeroids[ki], 
                            cen_bg1_qua, cen_bg1_cov_low, cen_bg1_cov_high, centroids[ki][4], 
                            cen_bg2_qua, cen_bg2_cov_low, cen_bg2_cov_high, centroids[ki][7], 
                            cluster_counts[ki]);
                }
                else
                {
                    sprintf(real_value_centeroids[ki], "    cluster %2d:", ki);
                    for(int noni=0; noni<m; noni++)
                    {
                        sprintf(real_value_centeroids[ki], "%s ~~~", real_value_centeroids[ki]);
                    }
                    sprintf(real_value_centeroids[ki], "%s %ld\0", 
                            real_value_centeroids[ki], 
                            cluster_counts[ki]);
                }
                printf("\t%s \n", real_value_centeroids[ki]);
            } 
        }
        else if(strcatCMD.find("SHOREmap backcross")!=std::string::npos && m>=2 && m<=5)
        {
            if(verbose) printf("Non-scaled centers of clusters given in PDF.\n");
            double cen_fg_cov_low, cen_fg_cov_high, cen_fg_qua;
            double cen_bg_cov_low, cen_bg_cov_high, cen_bg_cov, cen_bg_qua;
            
            for(int ki=k-1; ki>=0; ki --)
            {
                cen_fg_cov_low  = (1-centroids[ki][0])*(double)cluster_dim_coverage;
                cen_fg_cov_low  = cen_fg_cov_low*2.0;
                cen_fg_cov_high = (double)cluster_avg_coverage + cen_fg_cov_low;
                cen_fg_cov_low  = (double)cluster_avg_coverage - cen_fg_cov_low;
                if(cen_fg_cov_low < (double)filter_min_coverage)
                {
                    cen_fg_cov_low = cen_fg_cov_high;
                }
                cen_fg_qua      = centroids[ki][1]*max_quality;
            
                if(m == 4)
                {
                    cen_bg_cov      = (1-centroids[ki][2])*(double)bg_read;
                    cen_bg_qua      = max_quality-centroids[ki][3]*max_quality;
                }
                else if(m == 5)
                {
                    cen_bg_cov_low  = (1-centroids[ki][2])*bgref_cdim;
                    cen_bg_cov_low  = cen_bg_cov_low*2.0;
                    cen_bg_cov_high = bgref_cmean + cen_bg_cov_low;              
                    cen_bg_cov_low  = bgref_cmean - cen_bg_cov_low;
                    if((int)round(cen_bg_cov_low) < bg_ref_cov)
                    {
                        cen_bg_cov_low = cen_bg_cov_high;
                    }
                    cen_bg_qua      = centroids[ki][3]*max_quality;
                }
                else ;
                
            
                if(cluster_counts[ki] > 0)
                {
                    sprintf(real_value_centeroids[ki], "    cluster %2d: (L:%3.1f, H:%3.1f) %3.1f", ki,
                            cen_fg_cov_low, cen_fg_cov_high, cen_fg_qua);
                    if(m == 5)
                    {
                        sprintf(real_value_centeroids[ki],"%s (L:%3.1f, H:%3.1f) %3.1f %3.2f %8ld\0", 
                                real_value_centeroids[ki], 
                                cen_bg_cov_low, cen_bg_cov_high, cen_bg_qua, centroids[ki][4], 
                                cluster_counts[ki]);
                    }
                    else if(m ==4 )
                    {
                        sprintf(real_value_centeroids[ki], "%s %3.1f %3.1f %8ld\0", 
                                real_value_centeroids[ki],
                                cen_bg_cov, cen_bg_qua,
                                cluster_counts[ki]);
                    }
                    else if(m == 2)
                    {
                        sprintf(real_value_centeroids[ki], "%s %8ld\0", 
                                real_value_centeroids[ki],
                                cluster_counts[ki]);
                    };
                }
                else
                {
                    sprintf(real_value_centeroids[ki], "    cluster %2d:", ki);
                    for(int noni=0; noni<m; noni++)
                    {
                        sprintf(real_value_centeroids[ki], "%s ~~~", real_value_centeroids[ki]);
                    }
                    sprintf(real_value_centeroids[ki], "%s %ld\0", 
                            real_value_centeroids[ki], 
                            cluster_counts[ki]);
                }
		////printf("\t%s \n", real_value_centeroids[ki]);
            }
        }
        else
        {
            ;
        }
    }    
    /* figure out range of chromosome to plot ....................................................*/
    double ci_start = 0.0;
    double ci_end   = 0.0;
    double ci_step  = 0.0;
    if(strcatCMD.find("--chromosome") != std::string::npos)
    {
        ci_start = reg_begin;
        ci_end   = reg_end;
        ci_step  = (ci_end-ci_start+1)/5;
    }
    else
    {
        ci_start = 0.0; 
        ci_end   = (double)(CHR2SIZE[chrID]);
        //ci_step  = (double)ci_end/5;
        //(double)chrsizes_max
        //double ratio = (double)ci_end/(double)chrsizes_max;
        ci_step  = 3000000;
        if((double)chrsizes_max/ci_step < 1)
        while((double)chrsizes_max/ci_step < 5)                 // to avoid dense plotting of x-axis
        {
            ci_step -= 1000;
        }
        else
        while((double)chrsizes_max/ci_step > 10)                // to avoid dense plotting of x-axis
        {
            ci_step += 1000000;
        }
    }
    
    /* prepare window-averaged allele frequency, boost value, coefficient of variation........... */
    double winbeg;                        // position           of            window           begin
    double winend;                        // position           of            window             end
    double winctr;                        // position           of            window          center
    double avgfrq;                        // averaged                         window     frequency 1
    double avgbst;                        // boost value        of            window
    double bstmax;                        // maximum            of            boost           values
    double mut_sum;                       // coverage           of all mut-alleles     in  a  window
    double cov_sum;                       // coverage           of all ref+mut alleles in  a  window
    double sgl_frq;                       // allele frequency   of a      single                 SNP
    double sgl_sco;                       // base quality score of a      single     mutant     base
    long   mkr_sum;                       // number             of SNPs   in         a        window
    bstmax = 0.0;
    std::vector<WININFO> bstmax_pos;      // positions          of boost                         max
    std::vector<WININFO> winAVGnew;       // sum of counts(allele1) / sum of counts(allele1+allele2)
    std::vector<WININFO> winAF2new;       // sum of AFs(allele1+allele2) /num of markers in a window
    std::vector<WININFO> winSTDnew;       // std of  AFs   of     the     above          calculation         
    std::vector<WININFO> winBSTnew;       // boost value calculated from a  window-freq  of  markers
    std::vector<WININFO> winSCOnew;       // score value calculated for  a window-freq   of  markers    
    std::vector<double>  winSGLnew;       // AFs    of    single     markers     in     a     window
    
    if(plot_window)
    for(winbeg=(double)ci_start; winbeg<=(double)(ci_end-window_size+1); winbeg+=(double)window_step)
    {
        winend  = winbeg + window_size - 1;
        mut_sum = 0;
        cov_sum = 0;
        mkr_sum = 0;
        sgl_frq = 0;  
        sgl_sco = 0;
        winSGLnew.clear();                // clear tmp winSGLnew for next window....................
        mkr_itr     = interData.begin();
        mkr_itr_end = interData.end();
        while(mkr_itr != mkr_itr_end)
        {
            double mkrp = (double)(*mkr_itr).first;
            if(mkrp<winbeg)
            {
                mkr_itr++;
                continue;
            }
            if(mkrp>=winbeg && mkrp<=winend)
            {
                std::stringstream posTmp;
                posTmp << (unsigned long)(*mkr_itr).first;
                string qkey("");
                qkey += chrID;
                qkey += ".#.";
                qkey += posTmp.str();
                map<std::string, std::string>::iterator qitr = QUALITY1.find(qkey); 
                string qua = (*qitr).second.c_str();                  // like qua = eop002,40.......
                std::vector<string> quainfo = split_string(qua, ',');
                
                if((double)atol(quainfo[1].c_str()) < quality_min && 
                     strcatCMD.find("SHOREmap backcross")!=std::string::npos) 
                {
                    mkr_itr ++;
                    continue;
                }
                sgl_sco      += (double)(atol(quainfo[1].c_str()));
                
                TRIPLE cntTmp = (*mkr_itr).second;
                mut_sum      += (double)(cntTmp.Ci[1]);
                cov_sum      += (double)(cntTmp.Ci[0] + cntTmp.Ci[1] + cntTmp.Ci[2]);//caution 2014-05-05: with error
                mkr_sum      += 1;
                double tmp_sglfeq;
                tmp_sglfeq    = (double)(cntTmp.Ci[1])/(double)(cntTmp.Ci[0]+ cntTmp.Ci[1]+ cntTmp.Ci[2]);
                sgl_frq      += tmp_sglfeq;
                winSGLnew.push_back(tmp_sglfeq);
                
                mkr_itr ++;
            }
            if( mkr_itr==mkr_itr_end || mkrp>winend )
            {
                if(mkr_sum>=filter_min_marker)
                {
                    winctr   = (winend+winbeg)/2;
                    avgfrq   = mut_sum/cov_sum;
                    WININFO iwin;
                    iwin.pos = winctr;
                    iwin.beg = winbeg;                 //// is starting from the pos of the 1st snp?
                    iwin.end = winend;
                    iwin.val = avgfrq;
                    /* window average 1  */
                    winAVGnew.push_back(iwin);
                    /* std for winAF2new */
                    std::vector<double>::iterator sglfrq_itr;
                    std::vector<double>::iterator sglfrq_itr_end;
                    sglfrq_itr     = winSGLnew.begin();
                    sglfrq_itr_end = winSGLnew.end();
                    double x_xi_sq = 0;
                    while(sglfrq_itr != sglfrq_itr_end)
                    {
                        x_xi_sq += pow(*sglfrq_itr - iwin.val, 2);
                        sglfrq_itr ++;
                    }
                    iwin.val = sqrt(x_xi_sq/(winSGLnew.size()-1));
                    winSTDnew.push_back(iwin);
                    /* mean == winAF2new */
                    iwin.val = sgl_sco/(double)mkr_sum;
                    winSCOnew.push_back(iwin);
                    iwin.val = sgl_frq/(double)mkr_sum;
                    winAF2new.push_back(iwin);
                    
                    /*   boost  */
                    if(plot_boost)
                    {
                        avgfrq = sgl_frq/(double)mkr_sum;
                        avgbst = fabs(1 - max(expect, 1-expect)/max(avgfrq, 1-avgfrq));
                        if(avgbst == 0)
                        {
                            avgbst = (double)INF;
                        }
                        else
                        {
                            avgbst = 1/avgbst;
                        }

                        if(avgbst!=INF && avgbst>bstmax) 
                        {    
                            bstmax = avgbst;
                            if(bstmax_pos.size() > 0) bstmax_pos.clear();
                            iwin.val = bstmax;
                            bstmax_pos.push_back(iwin);
                        }
                        else if(avgbst==bstmax) 
                        {
                            iwin.val = bstmax;
                            bstmax_pos.push_back(iwin);
                        }
                        else ;
                        iwin.val = avgbst;
                        winBSTnew.push_back(iwin);
                    }
                }
                break;
            }
        }
    }
    /* check if windows exist to plot ............................................................*/
    //if(winAF2new.size() == 0   && strcatCMD.find("SHOREmap outcross")!=std::string::npos)
    if(plot_window && winAF2new.size() == 0)
    {
        printf("Warning: window-plot is asked, but it is cancelled (no boost-plot neither). \n");
        printf("Reason: there is not any window with min-window-markers=%ld. \n", 
                filter_min_marker);
        printf("Hint: increase window size or decrease min-window-markers from the cmd line.\n");
        plot_window  = false;
    }
    else
    if(plot_window && verbose)
    {
        printf("\tregion to plot: (chr_start, chr_end) = %.0f, %.0f \n", ci_start+1, ci_end);
        printf("\tnumber of windows-of-markers: %ld, ", winAF2new.size());
        printf("1st-center: %.0f, ", (*winAF2new.begin()).pos);
        std::vector<WININFO>::iterator avg_itr = --winAF2new.end();
        printf("last-center: %.0f\n", (*avg_itr).pos);
    }
    ////
    ////
    /* plot                             window-averaged                                 frequency */
    unsigned long  awinSize;
    double* winPosiSet = NULL;
    double* winFreqSet = NULL;
    double* winAFr2Set = NULL;
    double* winABSTSet = NULL;
    double* winSTDnew2 = NULL;
    double* winScorSet = NULL;
    if(plot_window)
    {
        /* prepare.........double*.........variables.........required............by........dislin */
        awinSize = winAVGnew.size();
        winPosiSet = (double*)malloc((awinSize+1)*sizeof(double));
        winFreqSet = (double*)malloc((awinSize+1)*sizeof(double));
        winAFr2Set = (double*)malloc((awinSize+1)*sizeof(double));
        winSTDnew2 = (double*)malloc((awinSize+1)*sizeof(double));
        winScorSet = (double*)malloc((awinSize+1)*sizeof(double));        
        if(plot_boost)
        {
            winABSTSet = (double*)malloc((awinSize+1)*sizeof(double));
        }
        if(winPosiSet==NULL || winFreqSet==NULL || (plot_boost && winABSTSet==NULL) || 
           winScorSet==NULL)
        {
            printf("Malloc error 2 in plot_chr_winboost(...). Exited.\n");
            exit(1);
        }
       /* prepare........xray:position.......and.......yray:frequency........for.........plotting */
       /* window-averaged.............................................................frequency 1 */
       std::vector<WININFO>::iterator win_itr     = winAVGnew.begin();
       std::vector<WININFO>::iterator win_itr_end = winAVGnew.end();  
       unsigned long                  winpos      = 0;
       while(win_itr != win_itr_end)
       {
            /* this position */
            *(winPosiSet+winpos) = (*win_itr).pos;
            *(winFreqSet+winpos) = (*win_itr).val;
            /* next position */
            win_itr ++;
            winpos ++;
        }
        /* std.........related............to.............window-averaged..............frequency 2 */
        win_itr     = winSTDnew.begin();
        win_itr_end = winSTDnew.end();
        winpos      = 0;
        while(win_itr != win_itr_end)
        {
            winSTDnew2[winpos] = (*win_itr).val;
            win_itr ++;
            winpos ++;
        }
        /* window-averaged............................................................frequency 2 */
        win_itr     = winAF2new.begin();
        win_itr_end = winAF2new.end();
        winpos      = 0;
        while(win_itr != win_itr_end)
        {
            *(winAFr2Set+winpos) = (*win_itr).val;
            win_itr ++;
            winpos  ++;
        }
        // winABSTSet //
        /* window-averaged..................................................................score */
        win_itr     = winSCOnew.begin();
        win_itr_end = winSCOnew.end();
        winpos      = 0;
        while(win_itr != win_itr_end)
        {
            *(winScorSet+winpos) = (*win_itr).val;
            win_itr ++;
            winpos  ++;
        }
        /* finding....a....region....around.......the.......lowest.......cv=std/mean.......values */
        if(strcatCMD.find("SHOREmap outcross")!=std::string::npos)
        {
            interest_regPos[0] = INF;
            interest_regPos[1] = 0;
            interest_regHET[0] = 0.025;
            interest_regHET[1] = 0.025;
            win_itr     = winSTDnew.begin();
            win_itr_end = winSTDnew.end();
            unsigned long gi = 0;
            while(win_itr != win_itr_end)
            {
                if((*win_itr).val <= cv_threshold)
                if(winAFr2Set[gi] >= avg_threshold && winAFr2Set[gi] <= interval_max_mean)
                {
                    /* record....minimum.....window....begin....and....maximum....window......end */
                    if((*win_itr).beg < interest_regPos[0]) interest_regPos[0] = (*win_itr).beg;
                    if((*win_itr).end > interest_regPos[1]) interest_regPos[1] = (*win_itr).end;
                }
                gi ++;
                win_itr ++;
            }
            /* if.........there............is..............a.............valid.............region */
            if(interest_regPos[1] > interest_regPos[0])
            {
                printf("\tMAPPING INTERVAL PREDICTED FROM %.0f ", interest_regPos[0]);
                printf("TO %.0f, ", interest_regPos[1]);
                if((unsigned long)(interest_regPos[1]-interest_regPos[0]+1)/1000000 > 0)
                {
                    sprintf(iregion, "*Predicted mapping interval of size %.2f Mbp (large): %.0f ~ %.0f\0", 
                            (interest_regPos[1]-interest_regPos[0]+1)/1000000.0,
                            interest_regPos[0], interest_regPos[1]);
                    printf("REGION SIZE = %.0f Mbp.\n",(interest_regPos[1]-interest_regPos[0]+1)/1000000.0);
                }
                else
                {
                    sprintf(iregion, "*Predicted mapping interval of size %.2f Kbp (normal): %.0f ~ %.0f\0", 
                            (interest_regPos[1]-interest_regPos[0]+1)/1000.0,
                            interest_regPos[0], interest_regPos[1]);
                    printf("REGION SIZE = %.0f Kbp.\n",(interest_regPos[1]-interest_regPos[0]+1)/1000.0);
                }
            }
        }
    }
    
    /* replace    boost    values     of     INF     as     the      maixmum     boost     value  */
    if(plot_window && plot_boost && strcatCMD.find("SHOREmap outcross")!=std::string::npos) 
    {
        std::vector<WININFO>::iterator bst_itr;
        std::vector<WININFO>::iterator bst_itr_end;
        bst_itr     = winBSTnew.begin();
        bst_itr_end = winBSTnew.end();
        while(bst_itr != bst_itr_end)
        {
            if((*bst_itr).val == INF)
            {
                (*bst_itr).val = bstmax;
                WININFO iwin;
                iwin.pos = (*bst_itr).pos;
                iwin.val = bstmax;
                bstmax_pos.push_back(iwin);
            }
            (*bst_itr).val = (*bst_itr).val/(bstmax+MARGIN);// normalized.........as...........[0,1]
            bst_itr ++;
        }
        
        /* new                              on                                  2013-07-16 18:22  */
        /* if the minimum boost value  is  around  0.5,  then  boost  value  is  not  informative */
        /* then we can pull the values down a little bit according to  the  minimum  boost  value */
        
        double bst_min = 1.0;
        bst_itr        = winBSTnew.begin();
        while(bst_itr != bst_itr_end)
        {
            if(bst_min  > (*bst_itr).val)
            {
                bst_min = (*bst_itr).val;
            }
            bst_itr ++;
        }
        if(bst_min >= 0.25 && interest_regPos[1] < interest_regPos[0])
        {
            // CAUTION: only pull down boost-values of a chromosome  indicating  no  confident  peak
            // TODO: also  consider   mean   of    AFs;    only    when    AF_max    <    say    0.6 
            bst_itr     = winBSTnew.begin();
            while(bst_itr != bst_itr_end)
            {
                (*bst_itr).val -= bst_min;
                bst_itr ++;
            }
        }
        /* window-boost.....................................................................value */
        if(plot_boost)
        {
            std::vector<WININFO>::iterator win_itr;
            std::vector<WININFO>::iterator win_itr_end;  
            win_itr     = winBSTnew.begin();
            win_itr_end = winBSTnew.end();
            unsigned long winpos = 0;
            while(win_itr != win_itr_end)
            {
                *(winABSTSet+winpos) = (*win_itr).val;        
                win_itr ++;
                winpos  ++;
            }
        }
    }
    
    /* output.......statistical.......info.......if.........required.........--........2013-12-23 */
    if(plot_record)
    {
        print_plot_info(chrID, myPosiSet,  myFreqSet, dtSize,
                               winPosiSet, winAFr2Set, winSTDnew2, winABSTSet, 
                               awinSize);
    }
    /* setting of  page  format, file   forma   and    file    name   in    the   parent-function */              
    /* set                                    axis                                         system */
    double unit_color;
    double base_color;
    double peak_color;
    if(strcatCMD.find("SHOREmap outcross")!=std::string::npos && k >= 2)
    {
        ////base_color = 254.0*0.35;
        ////peak_color = 254.0*0.55;
        base_color = 1.0;
        peak_color = 254.0;
        unit_color = (peak_color-base_color)/(k-1);
    }
    else if(k >= 2)
    {
        base_color = 1.0;
        peak_color = 254.0;
        unit_color = (peak_color-base_color)/(k-1);
    }     
    else ;       
    axspos(700,3500);                      // level 1    - determines position  of  an  axis  system
    int xaxisLen = 9700;
    int yaxisLen = 2200;
    int zaxisLen = 2200;
    if(strcatCMD.find("--chromosome") != std::string::npos)
    {
        ax3len(xaxisLen, yaxisLen, zaxisLen);
        if(k >= 2) colran (base_color, peak_color);
    }
    else                                 // scale  with   the   maximum   length   of    chromosomes
    {
        if(plot_scale)
        {
            double ratio = (double)ci_end/(double)chrsizes_max;
            ////axslen((int)(10100.0*ratio),2200);
            ax3len((int)(xaxisLen*ratio), yaxisLen, zaxisLen);
            if(k >= 2) colran (base_color, peak_color);
        }
        else
        {
            ////axslen(10100,2200);
            ax3len(xaxisLen, yaxisLen, zaxisLen);
            if(k >= 2) colran (base_color, peak_color);
        }
    }
    shdmod ("SYMB", "CURVE");
    //setscl(myPosiSet, dtSize, "x");      // level 1     - sets     the     scale      of      axis
    //setscl(myFreqSet, dtSize, "y");      // level 1    
    pagera();                              // level 1/2/3 - plot   a   border   around   the    page
    complx();                              // level       - complex                             font                             
    int ic0 = intrgb(0,0,0);               // level 1/2/3 - creates explicit color  value  from  RGB
    frmclr(ic0);                           // level 1/2/3 - defines      color       of       frames
    axclrs(ic0, "ALL", "XYZ");             //               ’LINE’, ’TICKS’, ’LABELS’, ’NAME’, ’ALL’
    height(80);                            // level 1/2/3 - defines height  of  characters  in  plot 
                                           //               (names  of  title&axis   not   included)
    helve();
    psfont("Helvetica");                                           
    name("Chromosome Position", "x");      // level 1/2/3 - defines           axis            titles
    if(plot_boost && plot_window)          //
    name("Allele Frequency (with boost)", "y");
    else
    name("Allele Frequency", "y");
    if(k >= 2) 
    name("Marker cluster", "z");
    hname(80);                             // level 1/2/3 - defines character height for axis  names
    labdig(-1, "x");                       // level 1/2/3 - defines number of decimal plcs in labels
    ticks(5, "x");                         // level 1/2/3 - defines number of ticks  between  labels
    ticks(1, "y");                         //
    /* set                                                                                  title */
    std::string myTitle = "Chromosome " + chrID;
    std::stringstream chrlength;
    chrlength << (unsigned long)CHR2SIZE[chrID];
    myTitle += ": ";
    myTitle += chrlength.str();
    myTitle += "bp";
    titlin(myTitle.c_str(), 1);            // level 1/2/3 - defines up to four lines  of  text  used 
                                           //               for       axis       system       titles
    htitle(80);                            // level 1/2/3 - defines  character  height  for   titles
                                           //               The character height defined  by  HEIGHT 
                                           //               will be used if  HTITLE  is  not  called  
    if(plot_window && plot_boost) 
    {
        if(verbose) printf("\tbstmax = %.3f @position(s): ", bstmax);
        
        std::string peak_file = "";
        peak_file = out_folder+ "SHOREmap_boosted_peaks.txt\0"; 
        fstream fpeakout(peak_file.c_str(), ios::out | ios::app);
        if(!fpeakout.is_open())
        {
            printf("Cannot open file to write boost peaks (in plot_chr_winboost.cpp). ");
            printf("Exited. \n");
            exit(1);
        }
        std::vector<WININFO>::iterator bst_max_itr;
        std::vector<WININFO>::iterator bst_max_itr_end;        
        bst_max_itr     = bstmax_pos.begin();
        bst_max_itr_end = bstmax_pos.end();
        int bsto = 1;
        while(bst_max_itr != bst_max_itr_end)
        {
            if(verbose && bsto <= 10)
            {
                printf("%ld \t", (unsigned long)(*bst_max_itr).pos);
                bsto ++;
            }
            fpeakout << chrID <<  "\t" << (unsigned long)(*bst_max_itr).pos << endl; // boosted-peak
            bst_max_itr ++;
        }
        if(verbose) printf("\n");
        if(verbose && bsto==10 && bstmax_pos.size()>10)
        {
            printf("%ld more positions with bstmax but not shown. \n", bstmax_pos.size()-10);
        }
        fpeakout.close();
    }
    
    if(k >= 2)  // 3D-plot..........................................................................
    {                             
        double rank_step = 1.0;
        //if(clusterK >= 20)
        if(k >= 20)
        {
            rank_step = 4.0;
        }
        colran (base_color, peak_color);
        labdig (-1, "Z");
        ////axsbgd (0);
        graf3(ci_start,      ci_end, ci_start,    ci_step,
                   0.0,        1.05,      0.0,        0.1,
                     0,         k-1,        0,  rank_step);
    }
    else      // 2D-plot....................................................................outcross
    {
        graf(ci_start,      ci_end, ci_start, ci_step,
                  0.0,        1.05,      0.0,    0.1);   // level 1 - plots..........2D.........axis
                                                         // system: x,y-lower/upper, 1-st-label/step
    }                       
    
    /* plot   allele-frequency    with    indications    from    parents    -    2013-06-14 18:41 */
    bool bgs_separated = false;
    for(unsigned long afi = 0; afi < dtSize; afi ++)
    {
        std::stringstream ss;
        ss << (unsigned long)myPosiSet[afi];
        std::string fndkey = (chrID + ".#." + ss.str());
        if(QUALITY1[fndkey].find("FLAG4parentB") != std::string::npos)
        {
            /* plot....................two.................background..................separately */
            bgs_separated = true;
            break;
        }
    }
    bgs_separated = false; // caution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! turned off on 2014-03-28 21:54
    if(bgs_separated && strcatCMD.find("SHOREmap outcross")!=std::string::npos)
    {
        /* check the number of SNPs from background parent a and b............................... */
        unsigned long SizePA = 0;
        unsigned long SizePB = 0;
        for(unsigned long afi = 0; afi < dtSize; afi ++)
        {
            std::stringstream ss;
            ss << (unsigned long)myPosiSet[afi];
            std::string fndkey = (chrID + ".#." + ss.str());
            if(QUALITY1[fndkey].find("FLAG4parentB") != std::string::npos)
            {
                SizePB ++;
            }
            else
            {
                SizePA ++;
            }
        } 
        /* prepare double array for dislin .......................................................*/
        double* myPosiSetPA = (double*)malloc((SizePA+1)*sizeof(double));
        double* myFreqSetPA = (double*)malloc((SizePA+1)*sizeof(double));
        double* myScorSetPA = (double*)malloc((SizePA+1)*sizeof(double)); 
        double* myPosiSetPB = (double*)malloc((SizePB+1)*sizeof(double));
        double* myFreqSetPB = (double*)malloc((SizePB+1)*sizeof(double));
        double* myScorSetPB = (double*)malloc((SizePB+1)*sizeof(double));         
        if(myPosiSetPA==NULL || myFreqSetPA==NULL || myPosiSetPB==NULL || myFreqSetPB==NULL ||
           myScorSetPA==NULL || myScorSetPB==NULL)
        {
            printf("Malloc error in plot_chr_winboost(...). Exited.\n");
            exit(1);
        }
        unsigned long iPA = 0;
        unsigned long iPB = 0;
        for(unsigned long afi = 0; afi < dtSize; afi ++)
        {
            std::stringstream ss;
            ss << (unsigned long)myPosiSet[afi];
            std::string fndkey = (chrID + ".#." + ss.str());
            if(QUALITY1[fndkey].find("FLAG4parentB") != std::string::npos)
            {
                myPosiSetPB[iPB] = myPosiSet[afi];
                myFreqSetPB[iPB] = myFreqSet[afi];  
                myScorSetPB[iPB] = myScorSet[afi];              
                iPB ++;
            }
            else
            {
                myPosiSetPA[iPA] = myPosiSet[afi];
                myFreqSetPA[iPA] = myFreqSet[afi];
                myScorSetPA[iPA] = myScorSet[afi];                              
                iPA ++;
            }

        } 
        /* plot       the                double      array              for                dislin */
        /* from       background         parent      b         -        without         phynotype */
        // color("magenta");
        solid();                                     //             - sets  a   solid   line   style
        hsymbl(25);                                  // level 1/2/3 - defines   size   of    symbols 
        // color("cyan");
        int ic;
        ic = intrgb(0.6835938,0.8398438, 0.9257812); // level 1/2/3 - creates explicit  color  value
        setclr(ic);        
        incmrk(-1);                                  // level 1/2/3 - selects symbol mode for  CURVE              
        marker(18); 
        
        if(plot_marker)
        if(k>=2)
        {
            curve3(myPosiSetPB, myFreqSetPB, myScorSetPB, SizePB);
        }
        else
        {
            curve(myPosiSetPB, myFreqSetPB, SizePB);
        }
        /* from       background         parent       a          -          with       phynotype  */
        solid();
        hsymbl(25); 
        //color("gray");
        ic = intrgb(0.6835938,0.8398438, 0.9257812);
        setclr(ic);        
        incmrk(-1);     
        marker(15);
        if(plot_marker)
        if(k>=2)
        {
            curve3(myPosiSetPA, myFreqSetPA, myScorSetPA, SizePA);   
        } 
        else
        {
            curve(myPosiSetPA, myFreqSetPA, SizePA);
        }  
        free(myPosiSetPA);
        free(myFreqSetPA);
        free(myScorSetPA);
        free(myPosiSetPB);
        free(myFreqSetPB);  
        free(myScorSetPB);
    }
    else
    {
        // PCI plot for backcross only
        if(pci && pci_chr.compare(chrID)==0)
        {
             solid();
             linwid(0.5);
             color("blue");
             marker(15);
             errbar(PciPosiSet, PciFreqSet, PciElowSet, PciEuppSet, ipci);
             incmrk(0);
             double xbe[2];
             double ybe[2];
             xbe[0] = pci_start;
             xbe[1] = pci_end;
             ybe[0] = oberveredAFMean;
             ybe[1] = oberveredAFMean;
             color("red");
             curve(xbe, ybe, 2);             
        }
        if(verbose && plot_marker)
        printf("\tNumber of markers to plot: %ld\n", dtSize);
        if(strcatCMD.find("SHOREmap outcross")!=std::string::npos)
        {
            hsymbl(40);
        }
        else
        {
            hsymbl(25);
        }
        marker(15);
        if(plot_marker)
        if(k >= 2) 
        {
            curve3(myPosiSet, myFreqSet, myScorSet, (int)dtSize);   
        }
        else
        {
            solid();
            incmrk(-1);
            color("red");
            //int ic = intrgb(0.6835938,0.8398438, 0.9257812);
            //setclr(ic);
            curve(myPosiSet, myFreqSet, dtSize);         // level 0/1   - marks points with symbols.
        }
    }
    /* used before scatter plot                                                                   */
    // curve(myPosiSet, myFreqSet, dtSize);// level 0/1   - marks data points with symbols.
    /* end of - 2013-06-14 18:41                                                                  */
    /* set markers                         */
    marker(1);                             // level 1/2/3 - select symbols  used  to   plot   points 
                                           //               The  symbol number  will  be incremented 
                                           //               by 1 after a certain number of calls  to
                                           //               CURVE       defined      by       INCCRV
    hsymbl(30);                            // level 1/2/3 - defines      size       of       symbols 
    /* background                                                                           color */
    // int ic = intrgb(0.85,0.85,0.85);    // level 1/2/3 - create  explicit color  value  from  RGB
    // axsbgd (ic);                        // level 1/2/3 - define  background color of axis systems
    setrgb(0.7, 0.7, 0.7);                 // level 1/2/3 - define   foreground   color   from   RGB
    /* plotting...............................data........................................points  */
    if(plot_window && strcatCMD.find("SHOREmap backcross")!=std::string::npos)
    {
        if(winAVGnew.size()==1)
        {
            incmrk(-1);
        }
        else 
        {
            incmrk(100); 
            dashl();
        }
        /* window-averaged .................. frequency........................................AF2*/             
        marker(3); 
        int ic = intrgb(0, 0, 1);
        setclr(ic);
        curve(winPosiSet, winAFr2Set, awinSize);
    }
    if(plot_window && strcatCMD.find("SHOREmap outcross")!=std::string::npos)
    {
        /* window-averaged frequency - AF 1 -- NOT plotted after 2013-10-15 15:07
            marker(2);
            int ic = intrgb(0.9960938, 0.2578125, 0.3632812);
            setclr(ic);
            if(winAVGnew.size()==1) incmrk(-1);
            else incmrk(60);
            curve(winPosiSet, winFreqSet, awinSize);
        */   
        if(winAVGnew.size() == 1) 
        {
            incmrk(-1);
        }
        else 
        {
            incmrk(100); 
            dashl();
        }  
        /* window-averaged...................frequency.......................................AF 2 */             
        marker(3); 
        int ic = intrgb(0, 0, 1);
        setclr(ic);
        thkcrv(3);
        curve(winPosiSet, winAFr2Set, awinSize); 
        /* boost values...........................................................................*/
        if(plot_boost)
        {
            incmrk(160);
            color("gray");
            curve(winPosiSet, winABSTSet, awinSize);
        }
        /* mapping interval.......................................................................*/
        if(interest_regPos[1] > interest_regPos[0])
        {
            hsymbl(50);
            solid();
            incmrk(0);
            ic = intrgb(0.359375, 0.6523438, 0.7265625);
            setclr(ic); 
            thkcrv(30);
            curve(interest_regPos, interest_regHET, 2);
        }
    }
    int ict = intrgb(0, 0, 0);
    setclr(ict);
    title();                               // level 2/3   -   plots  title over   an   axis   system 
                                           //                 The title may contain up to four lines 
                                           //                         of text designated with TITLIN 
                                           //                           not required for qplsca(...)
    int ic = intrgb(0.54, 0.54, 0.54);     // level 1/2/3 -   creates explicit color value from  RGB
    setclr(ic);                            //                
    dotl();                                //             -   sets    a      dotted     line   style
    grid(0,1);                             // level 2/3   -   plot..............................grid
    if(plot_window  && strcatCMD.find("SHOREmap outcross")!=std::string::npos)
    {
        ic = intrgb(0, 0, 0);
        setclr(ic);
        char windowINFO[512];
        if(CMD.find("--min-coverage")!=CMD.end() || CMD.find("--max-coverage")!=CMD.end())
        {
            std::string strparas(" ");
            strparas += "interval-mean-min~max \%.4f~\%.4f ";
            strparas += "with Cv-max \%.4f; ";
            strparas += "window-size \%ld ";
            strparas += "with step \%ld; ";
            strparas += "marker-min \%ld; ";
            if(CMD.find("--min-coverage")!=CMD.end() && CMD.find("--max-coverage")!=CMD.end())
            {
                strparas += "coverage-min~max \%ld~\%ld\0";
                sprintf( windowINFO, strparas.c_str(),
                         interval_min_mean, interval_max_mean, interval_max_cvar,
                         window_size, window_step, 
                         filter_min_marker, 
                         filter_min_coverage, filter_max_coverage); 
            }
            else
            if(CMD.find("--min-coverage") != CMD.end())
            {
                strparas += "coverage-min \%ld\0";
                sprintf( windowINFO, strparas.c_str(), 
                         interval_min_mean, interval_max_mean, interval_max_cvar,
                         window_size, window_step, 
                         filter_min_marker, 
                         filter_min_coverage);
            }
            else
            {
                strparas += "coverage-max \%ld\0";
                sprintf( windowINFO, strparas.c_str(), 
                         interval_min_mean, interval_max_mean, interval_max_cvar,
                         window_size, window_step, 
                         filter_min_marker, 
                         filter_max_coverage);
            }
        }
        else
        {
            std::string strparas(" ");
            strparas += "interval-mean-min~max \%.3f~\%.3f; ";
            strparas += "with Cv-max \%.3f; ";
            strparas += "window-size \%ld ";
            strparas += "with step \%ld; ";
            strparas += "marker-min \%ld\0";
            sprintf( windowINFO, strparas.c_str(), 
                     interval_min_mean, interval_max_mean, interval_max_cvar,
                     window_size, window_step, 
                     filter_min_marker);
        }
        if(plot_boost)
        {
            messag("*Parameter settings: ", 700, 4000);
            messag(windowINFO, 700, 4150);
        }
        else
        {
            messag("*Parameter settings: ", 700, 4000);
            messag(windowINFO, 700, 4150);
        }
    }
    /*...................................adding legends...........................................*/
    string tmp("");
    // PCI plot for backcross only
    if(strcatCMD.find("SHOREmap backcross")!=std::string::npos && pci && pci_chr.compare(chrID)==0)
    {
        tmp  += "*Red line: mean of AF within the given mapping interval; "; 
        tmp  += "*blue bar: Wilson score interval @wikipedia: binomial proportion confidence interval";
        messag(tmp.c_str(), 700, 4050);
        tmp.clear();
    }
    if(plot_marker)
    {
        if(k == 1)
        {
            tmp += "*Points in colors: AF (=alt divided by (alt+ref)) at markers. ";
        }
        else
        {
            tmp += "*Points in colors: AF (=alt divided by (alt+ref)) at markers with clustering. ";
        }
        tmp += "alt or ref: coverage of non-reference or reference allele";
        messag(tmp.c_str(), 700, 4350);
    }
    int ypos  = 4350;
    if(plot_window)
    {
        tmp.clear();
        tmp  += "*Dashed line in blue: window-based AF (= summation-of-single-marker-AF ";
        tmp  += "divided by number-of-markers with min quality score involved).";
        int ic = intrgb(0, 0, 1);
        setclr(ic);
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
        if(strcatCMD.find("SHOREmap outcross")!=std::string::npos)
        {
            if(plot_boost)
            {
                color("gray");
                ypos  += 150;
                messag("*Dashed line in gray: window-based boost value ", 700, ypos);
            }
            if(interest_regPos[1] > interest_regPos[0])
            {
                ic = intrgb(0.359375, 0.6523438, 0.7265625);
                setclr(ic);
                ypos  += 150;
                messag(iregion, 700, ypos);
            }
        }
    }
    if(k >= 2 && plot_marker)
    {
        int ic = intrgb(0, 0, 0);
        setclr(ic);
        if(m == 2)
        {
            tmp.clear();
            tmp   = "*rank/cluster table: min_fg_cov max_fg_cov fg_qua ";
            tmp  += "size_cluster";
            ypos += 150;
            messag(tmp.c_str(), 700, ypos);
        }
        else if(m == 4)
        {
            tmp.clear();
            tmp   = "*rank/cluster table: min_fg_cov max_fg_cov fg_qua ";
            tmp  += "bg_mut_cov bg_mut_qua ";
            tmp  += "size_cluster";
            ypos += 150;
            messag(tmp.c_str(), 700, ypos);
        }
        else if(m == 5)
        {
            tmp.clear();
            tmp   = "*rank/cluster table: min_fg_cov max_fg_cov fg_qua ";
            tmp  += "min_bg_ref_cov max_bg_ref_cov bg_ref_qua ";
            tmp  += "size_cluster";
            ypos += 150;
            messag(tmp.c_str(), 700, ypos);
        }
        else if(m==8)
        {
            tmp.clear();
            tmp  = "* rank/cluster table: (min_fg_cov, max_fg_cov) fg_qua ";
            tmp += "bg1_qua (min_bg1_cov, max_bg1_cov) bg1_af ";
            tmp += "bg2_qua (min_bg2_cov, max_bg2_cov) bg2_af ";
            tmp += "size_cluster";
            ypos += 150;
            messag(tmp.c_str(), 700, ypos);
        }
        else ;
        ypos += 250;
        int num_ypos = 0;
        for(int ranki = k-1; ranki >= 0; ranki --)
        {
            if(ranki > 0)
            {
                setclr(round(base_color + ranki*unit_color));
            }
            else
            {
                setclr(round(base_color));
            }
            if(ranki>k-1-20 || k<=20)
            {
                messag(real_value_centeroids[ranki], 700, ypos);
                num_ypos ++;
            }
            else
            {
                if(strcatCMD.find("SHOREmap outcross")!=std::string::npos)
                {
                    messag(real_value_centeroids[ranki], 5500, ypos-num_ypos*150);
                }
                else
                {
                    messag(real_value_centeroids[ranki], 4200, ypos-num_ypos*150);
                }
            }
            ypos += 150;
        }
    } 
    /* ...................................... termination ........................................*/
    errmod("PROTOCOL", "OFF");             // level 1/2/3 - disable  printing  protocol  -   caution
    endgrf();                              // level 2/3   - terminates an axis system; set  level  1
    /* release memory ............................................................................*/
    if(myPosiSet) free(myPosiSet);
    if(myFreqSet) free(myFreqSet);
    if(myScorSet) free(myScorSet);
    if(pci && pci_chr.compare(chrID)==0)
    {
        if(PciPosiSet) free(PciPosiSet);
        if(PciFreqSet) free(PciFreqSet);
        if(PciEuppSet) free(PciEuppSet);
        if(PciElowSet) free(PciElowSet);
    }
    if(plot_window  && strcatCMD.find("SHOREmap outcross")!=std::string::npos)
    {
        if(winPosiSet) free(winPosiSet);
        if(winFreqSet) free(winFreqSet);
        if(winAFr2Set) free(winAFr2Set);
        if(winSTDnew2) free(winSTDnew2);
        if(plot_boost)
        if(winABSTSet) free(winABSTSet);
    }
    if(k >= 2)
    {
        for(int mali=0; mali<dtSize+1; mali++)
        {
           if(markerData[mali]) free(markerData[mali]);
           if(mali <= k)
           if(centroids[mali])  free(centroids[mali]);
        }
        if(markerData) free(markerData);
        if(centroids)  free(centroids);
    
        for(int irvc=0; irvc<k+1; irvc++)
        {
            if(real_value_centeroids[irvc]) free(real_value_centeroids[irvc]);
        }
        if(real_value_centeroids) free(real_value_centeroids);
        if(cluster_counts) free(cluster_counts);
    }
    /* TODO - plot coverage on markers, 2013-09-17 12:52                                          */
    return true;
}
