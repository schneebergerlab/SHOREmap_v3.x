/* this function visualizes annotations of snps if asked by users */
#include                      <stddef.h>
#include                      <stdlib.h>
#include                       <stdio.h>
#include                        <math.h>
#include                        <string>
#include                      <string.h>
#include                           <map>
#include                       <sstream>
#include                      <iostream>
#include             "dislin/dislin_d.h"
#include                     "globals.h"
#include                 "read_marker.h"
#include          "read_allele_count2.h"
// string fchrsizes
// multimap<unsigned long, mySNP> SNPlist

bool plot_chr_snp_annotations()
{
    //  check specific mutations, get their number and set snps as markers to read consensus info
    multimap<unsigned long, mySNP>::iterator snpitr = SNPlist.begin();
    unsigned long num_CDSsyn     = 0;
    unsigned long num_CDSnonsyn  = 0;
    unsigned long num_sschange   = 0;
    unsigned long num_5pUTR      = 0;
    unsigned long num_3pUTR      = 0;
    unsigned long num_intergenic = 0;
    unsigned long num_intronic   = 0;
    unsigned long num_other      = 0;
    while(snpitr != SNPlist.end())
    {
        std::stringstream seqtype;
        seqtype.str("");
        seqtype << (*snpitr).second.stype;
        if((seqtype.str()).find("CDS") != std::string::npos)
        {
            if((*snpitr).second.ns_change == 1) // nonsyn
            {
                num_CDSnonsyn ++;
            }
            else                                // syn
            {
                num_CDSsyn ++;
            }
        }
        else if((seqtype.str()).find("splice_site_change") != std::string::npos)
        {
            num_sschange ++;
        }
        else if((seqtype.str()).find("five_prime_UTR") != std::string::npos)
        {
            num_5pUTR ++;
        }
        else if((seqtype.str()).find("three_prime_UTR") != std::string::npos)
        {
            num_3pUTR ++;
        }
        else if((seqtype.str()).find("intronic/noncoding") != std::string::npos)
        {
            num_intronic ++;
        }
        else if((seqtype.str()).find("intergenic") != std::string::npos)
        {
            num_intergenic ++;
        }
        else
        {
            num_other ++;
        }
        // set allele info for extracting consensus info
        std::stringstream ss;
        ss.str("");
        ss << reg_chromosome << ".#." << (*snpitr).second.position;
        ALLELE1.insert(std::pair<std::string, std::string>(ss.str(), (*snpitr).second.ref_base));
        ALLELE2.insert(std::pair<std::string, std::string>(ss.str(), (*snpitr).second.new_base));
        snpitr ++;
    }
    
    cout << "Summary of annotated mutations: " << endl;
    cout << "\tintergenic         mutations: "        << num_intergenic << endl;
    cout << "\tintronic/noncoding mutations: "        << num_intronic   << endl;
    cout << "\tthree_prime_UTR    mutations: "        << num_3pUTR      << endl;
    cout << "\tfive_prime_UTR     mutations: "        << num_5pUTR      << endl;
    cout << "\tsplice site change mutations: "        << num_sschange   << endl;
    cout << "\tCDS (synnonymous)  mutations: "        << num_CDSsyn     << endl;
    cout << "\tCDS (non-syn....)  mutations: "        << num_CDSnonsyn  << endl;
    cout << "\tother              mutations: "        << num_other      << endl;
    
    // read consensus base and allele1 allele2 error counts
    if (!read_allele_counts2((char*)fconsensus.c_str()))
    {
        printf("ERROR: no consensus info recorded. Exited. \n");
        cout << fconsensus << endl;
        exit(1);
        //   CHR2POS2_ale1_ale2_err_COUNT[chr][pos] : {count_allele2, count_allele1, count_error};
    }
    //
    if(CHR2POS2_ale1_ale2_err_COUNT.find(reg_chromosome) == CHR2POS2_ale1_ale2_err_COUNT.end())
    {
        printf("ERROR: no info of chromosome %s in the consensus file %s. Pls check your files. Visualization exited.", 
                reg_chromosome.c_str(), fconsensus.c_str());
        return false;
    }
    map<unsigned long, TRIPLE> internalData = CHR2POS2_ale1_ale2_err_COUNT[reg_chromosome];
    if(internalData.size() == 0) 
    {
        printf("No markers to visualize for chromosome %s.\n\n", reg_chromosome.c_str());
        return false;
    }
    
    // allocate variable for visualization
    double* Posi_CDSnonsyn  = (double*)malloc((num_CDSnonsyn+1)*sizeof(double));
    double* Freq_CDSnonsyn  = (double*)malloc((num_CDSnonsyn+1)*sizeof(double));
    double* Posi_CDSsyn     = (double*)malloc((num_CDSsyn+1)*sizeof(double));
    double* Freq_CDSsyn     = (double*)malloc((num_CDSsyn+1)*sizeof(double));    
    double* Posi_ss         = (double*)malloc((num_sschange+1)*sizeof(double));
    double* Freq_ss         = (double*)malloc((num_sschange+1)*sizeof(double));    
    double* Posi_5pUTR      = (double*)malloc((num_5pUTR+1)*sizeof(double));
    double* Freq_5pUTR      = (double*)malloc((num_5pUTR+1)*sizeof(double));     
    double* Posi_3pUTR      = (double*)malloc((num_3pUTR+1)*sizeof(double));
    double* Freq_3pUTR      = (double*)malloc((num_3pUTR+1)*sizeof(double));  
    double* Posi_intronic   = (double*)malloc((num_intronic+1)*sizeof(double));
    double* Freq_intronic   = (double*)malloc((num_intronic+1)*sizeof(double));    
    double* Posi_intergenic = (double*)malloc((num_intergenic+1)*sizeof(double));
    double* Freq_intergenic = (double*)malloc((num_intergenic+1)*sizeof(double)); 
    // get the (pos, allele-frequency) info for each snp mutation
    unsigned long i_CDSsyn     = 0;
    unsigned long i_CDSnonsyn  = 0;
    unsigned long i_sschange   = 0;
    unsigned long i_5pUTR      = 0;
    unsigned long i_3pUTR      = 0;
    unsigned long i_intergenic = 0;
    unsigned long i_intronic   = 0;
    snpitr = SNPlist.begin();
    while(snpitr != SNPlist.end())
    {
        std::stringstream seqtype;
        seqtype.str("");
        seqtype << (*snpitr).second.stype;
        unsigned long ipos = (*snpitr).second.position;
        map<unsigned long, TRIPLE>::iterator mkr_itr;
        mkr_itr = internalData.find(ipos);
        if(mkr_itr == internalData.end())
        {
            printf("WARNING: position %ld not found in the consensus file. Skipped.\n", ipos);
            snpitr ++;
            continue;
        }
        double iaf = (double)(*mkr_itr).second.Ci[0]+(*mkr_itr).second.Ci[1]+(*mkr_itr).second.Ci[2];
        iaf        = (double)(*mkr_itr).second.Ci[1]/iaf;           // ref/icov or mut/icov 
        
        if((seqtype.str()).find("CDS") != std::string::npos)
        {
            if((*snpitr).second.ns_change == 1) // nonsyn
            {
                *(Posi_CDSnonsyn + i_CDSnonsyn ) = ipos;
                *(Freq_CDSnonsyn + i_CDSnonsyn ) = iaf;
                i_CDSnonsyn ++;
            }
            else                                // syn
            {
                *(Posi_CDSsyn + i_CDSsyn ) = ipos;
                *(Freq_CDSsyn + i_CDSsyn ) = iaf;
                i_CDSsyn ++;
            }
        }
        else if((seqtype.str()).find("splice_site_change") != std::string::npos)
        {
            *(Posi_ss + i_sschange ) = ipos;
            *(Freq_ss + i_sschange ) = iaf;
            i_sschange ++;
        }
        else if((seqtype.str()).find("five_prime_UTR") != std::string::npos)
        {
            *(Posi_5pUTR + i_5pUTR ) = ipos;
            *(Freq_5pUTR + i_5pUTR ) = iaf;
            i_5pUTR ++;
        }
        else if((seqtype.str()).find("three_prime_UTR") != std::string::npos)
        {
            *(Posi_3pUTR + i_3pUTR ) = ipos;
            *(Freq_3pUTR + i_3pUTR ) = iaf;
            i_3pUTR ++;
        }
        else if((seqtype.str()).find("intronic/noncoding") != std::string::npos)
        {
            *(Posi_intronic + i_intronic ) = ipos;
            *(Freq_intronic + i_intronic ) = iaf;
            i_intronic ++;
        }
        else if((seqtype.str()).find("intergenic") != std::string::npos)
        {
            *(Posi_intergenic + i_intergenic ) = ipos;
            *(Freq_intergenic + i_intergenic ) = iaf;
            i_intergenic ++;
        }
        snpitr ++;
    }       
    
    char chrNameBuf[512];
    sprintf(chrNameBuf, "%ssnp_ann_visualization_chr%s_%ld_%ld.pdf\0", 
    out_folder.c_str(), reg_chromosome.c_str(), reg_begin, reg_end);
    metafl("pdf");                         // level 0     - defines    the      metafile      format
    setfil(chrNameBuf);                    // level 0     - sets         alternate          filename
    filmod("VERSION");                     // level 0/1/2/3 do not shorten existing filename
    scrmod("revers");                      // level 0     - swaps  back   and   foreground   colours
    setpag("USEL");                        // level 0     - selects   a   predefined   page   format
                                           //               landscape:       11180*8640       points
    /* initialization                                                                             */
    disini();                              // level 0     - initialize DISLIN with default PARAs and
                                           //               a plotfile. The level  is  set   to   1. 
                                           //               DISINI must be called before  any  other 
                                           //               DISLIN   routine     except    for   ...
    /* setting of  page  format, file   forma   and    file    name   in    the   parent-function */              
    /* set                                    axis                                         system */
    double unit_color;
    double base_color;
    double peak_color;                                         
    base_color = 1.0;
    peak_color = 254.0;
    unit_color = (peak_color-base_color)/7;            
    axspos(700,3500);                      // level 1    - determines position  of  an  axis  system
    int xaxisLen = 9700;
    int yaxisLen = 2200;
    int zaxisLen = 2200;
    //ax3len(xaxisLen, yaxisLen, zaxisLen);
    axslen(xaxisLen,yaxisLen);
    //shdmod ("SYMB", "CURVE");
    pagera();                              // level 1/2/3 - plot   a   border   around   the    page
    complx();                              // level       - complex                             font                             
    int ic0 = intrgb(0,0,0);               // level 1/2/3 - creates explicit color  value  from  RGB
    frmclr(ic0);                           // level 1/2/3 - defines      color       of       frames
    //axclrs(ic0, "ALL", "XYZ");           //               ’LINE’, ’TICKS’, ’LABELS’, ’NAME’, ’ALL’
    height(80);                            // level 1/2/3 - defines height  of  characters  in  plot 
                                           //               (names  of  title&axis   not   included)
    helve();
    psfont("Helvetica");                                           
    name("Chromosome Position", "x");      // level 1/2/3 - defines           axis            titles
    name("Allele Frequency",    "y");   
    hname(80);                             // level 1/2/3 - defines character height for axis  names
    labdig(-1, "x");                       // level 1/2/3 - defines number of decimal plcs in labels
    ticks(5, "x");                         // level 1/2/3 - defines number of ticks  between  labels
    ticks(1, "y");                         //
    /* set                                                                                  title */
    std::string myTitle = "Chromosome " + reg_chromosome;
    std::stringstream chrlength;
    chrlength.str("");
    chrlength << (unsigned long)CHR2SIZE[reg_chromosome];
    myTitle += ": ";
    myTitle += chrlength.str();
    myTitle += "bp";
    titlin(myTitle.c_str(), 1);            // level 1/2/3 - defines up to four lines  of  text  used 
                                           //               for       axis       system       titles
    htitle(100);                           // level 1/2/3 - defines  character  height  for   titles
                                           //               The character height defined  by  HEIGHT 
                                           //               will be used if  HTITLE  is  not  called
                                               
    double ci_start  = 0.0; 
    double ci_end    = (double)CHR2SIZE[reg_chromosome];
    double ci_step   = (ci_end-ci_start+1)/5;
    ci_start = reg_begin;
    ci_end   = reg_end;
    ci_step  = (ci_end-ci_start+1)/5;
    
    double rank_step = 1.0;
    graf(ci_start,   ci_end, ci_start, ci_step,
              0.0,     1.05,      0.0,    0.1);      // level 1 - plots..........2D.........axis
                                                     // system: x,y-lower/upper, 1-st-label/step
    color("gray");                         //                
    dotl();                                //             -   sets    a      dotted     line   style
    grid(0,1);                             // level 2/3   -   plot..............................grid
                                                         
    solid();                                         //             - sets  a   solid   line   style
    int ic;
    incmrk(-1);                                      // level 1/2/3 - selects symbol mode for  CURVE              
    string tmp("");
    int ypos  = 4350;
    std::stringstream ss;
    if(num_intergenic > 0)                    
    {
        marker(3);
        hsymbl(100);                                  // level 1/2/3 - defines   size   of    symbols                                                  
        color("GREEN");
        curve(Posi_intergenic, Freq_intergenic, (int)num_intergenic);
        tmp.clear();
        ss.str("");
        ss << num_intergenic;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in intergenic regions";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_intronic > 0)                    
    {
        marker(15);
        hsymbl(90);                                 // level 1/2/3 - defines   size   of    symbols                                                  
        color("GRAY");
        curve(Posi_intronic, Freq_intronic, (int)num_intronic);     
        tmp.clear();
        ss.str("");
        ss << num_intronic;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in intronic/noncoding regions";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_3pUTR > 0)
    {
        marker(19);
        hsymbl(100);                                  // level 1/2/3 - defines   size   of    symbols                                                          
        color("WHITE");
        curve(Posi_3pUTR, Freq_3pUTR, (int)num_3pUTR);  
        tmp.clear();
        ss.str("");
        ss << num_3pUTR;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in 3' UTR regions";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_5pUTR > 0)
    {
        marker(5);
        hsymbl(150);                                 // level 1/2/3 - defines   size   of    symbols                                                          
        color("BLUE");
        curve(Posi_5pUTR, Freq_5pUTR, (int)num_5pUTR);  
        tmp.clear();
        ss.str("");
        ss << num_5pUTR;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in 5' UTR regions";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_sschange > 0)
    {
        hsymbl(120);                                 // level 1/2/3 - defines   size   of    symbols                                                          
        color("MAGENTA");
        curve(Posi_ss, Freq_ss, (int)num_sschange);  
        tmp.clear();
        ss.str("");
        ss << num_sschange;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in splice sites";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_CDSsyn > 0)
    {
        hsymbl(90);                                 // level 1/2/3 - defines   size   of    symbols                                                          
        color("CYAN");
        curve(Posi_CDSsyn, Freq_CDSsyn, (int)num_CDSsyn);
        tmp.clear();
        ss.str("");
        ss << num_CDSsyn;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in CDS with synonymous effect";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    if(num_CDSnonsyn > 0)
    {
        marker(18);
        hsymbl(85);
        color("RED");
        curve(Posi_CDSnonsyn, Freq_CDSnonsyn, (int)num_CDSnonsyn);
        tmp.clear();
        ss.str("");
        ss << num_CDSnonsyn;
        tmp  += "*";
        tmp  += ss.str();
        tmp  += " mutations in CDS with nonsynonymous effect";
        ypos += 150;
        messag(tmp.c_str(), 700, ypos);
    }
    color("WHITE");
    tmp.clear();
    tmp  += "*LEGEND:";
    messag(tmp.c_str(), 700, ypos-1100);
    tmp.clear();
    tmp  += "*Note: multiple annotations for one mutation can exist due to ";
    tmp  += "multiple isoforms of one gene or overlaps of different genes.";    
    ypos += 300;
    messag(tmp.c_str(), 700, ypos);
    
    title();                               // level 2/3   -   plots  title over   an   axis   system     
    /* ...................................... termination ........................................*/
    errmod("PROTOCOL", "OFF");             // level 1/2/3 - disable  printing  protocol  -   caution
    endgrf();                              // level 2/3   - terminates an axis system; set  level  1                                       
    disfin();                              // level 1/2/3 - terminates DISLIN and prints a message on
                                           //               the screen. The level is set  back  to  0
                                           //               not    required      for      qplsca(...)             
    // release memory
    if(Posi_CDSnonsyn)  free(Posi_CDSnonsyn);
    if(Freq_CDSnonsyn)  free(Freq_CDSnonsyn);
    if(Posi_CDSsyn)     free(Posi_CDSsyn);
    if(Freq_CDSsyn)     free(Freq_CDSsyn);
    if(Posi_ss)         free(Posi_ss);
    if(Freq_ss)         free(Freq_ss);
    if(Posi_5pUTR)      free(Posi_5pUTR);
    if(Freq_5pUTR)      free(Freq_5pUTR);
    if(Posi_3pUTR)      free(Posi_3pUTR);
    if(Freq_3pUTR)      free(Freq_3pUTR);
    if(Posi_intronic)   free(Posi_intronic);
    if(Freq_intronic)   free(Freq_intronic);              
    if(Posi_intergenic) free(Posi_intergenic);
    if(Freq_intergenic) free(Freq_intergenic);  
    
    return true;
}
