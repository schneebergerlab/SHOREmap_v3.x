/* read gene feature information from GFF file                                                    */
/* note: from 2013-05-06 21:05 on, isoforms of a gene are handled as individual genes             */
/*       e.g., Aa_G9830.t1, Aa_G9830.t2, Aa_G9830.t3 are recorded separately in gene_ann          */

#include  <stdio.h>
#include <iostream>
#include <stdlib.h>
#include      <map>
#include <string.h>
#include <assert.h>
#include "globals.h"
#include "read_chromosome_seq.h"

bool read_GFFinfo(char* file_gff, std::string* isof)
{
    FILE* fp = fopen(file_gff, "r");
    if(!fp)
    {
        printf("Cannot open gff file. Exited.\n");
        exit(1);
    }
    /* Read reference chromosome sequence from file such as A_alpina_V3.fa.shore; caution below:  */
    /* reference file contains chr ids as "1","2","3",...; this can be case sensitive             */
    /* reg_chromosome.c_str()+6: "1", "2", "3", to be consistent with genome file                 */
    
    std::string RNAtype("mRNA"); // caution: there are ncRNA, snoRNA, tRNA etc.

    unsigned long max_chr_len;
    map<std::string, unsigned long>::iterator csize_itr;
    csize_itr   = CHR2SIZE.begin();
    max_chr_len = (*csize_itr).second;
    while(csize_itr != CHR2SIZE.end())
    {
        /* find maximum chromosome length */
        if((*csize_itr).second > max_chr_len)
        {
             max_chr_len = (*csize_itr).second;
        }
        csize_itr ++;
    }
    
    // caution: (char*)reg_chromosome.c_str()+6 is for alpina; because for alpina, chromosome ids in 
    // A_alpina_V3.chrSizes.txt and A_alpina_V3.annotation.gff3 are: Aa.chr1, Aa.chr2,...
    std::string chr_seq;
    std::string chridtemp = "";
    //if(reg_chromosome.length()>6) chridtemp = reg_chromosome.c_str()+6; // caution: TODO standard required
    //else chridtemp = reg_chromosome.c_str();
    
    chridtemp = reg_chromosome.c_str(); // 2013-08-02 09:24
    
    if(!read_chromosome_seq( (char*)genome_file.c_str(), 
                             (char*)chridtemp.c_str(), 
                             max_chr_len, 
                             &chr_seq) )
    {
        //printf("ERROR: cannot find chromosome %s ", (char*)reg_chromosome.c_str()+6);
        printf("ERROR: cannot find chromosome %s ", (char*)reg_chromosome.c_str());// 2013-08-02 09:36    
        printf("in %s. Exited (in function read_GFFinfo(...)). ", (char*)genome_file.c_str());
        printf("Hint: check if IDs of chromosomes in all files are consistent.\n");
        exit(1);
    }
    if(verbose) printf("Chromosome %s: length = %ld is ready.\n", 
        (char*)reg_chromosome.c_str(), chr_seq.length());
    
    unsigned long ssnum; // number of splice sites for each gene
    /* GFF format: 0.chr 1.src 2.seq_type 3.start 4.end 5.score 6.orien 7.frame 8.description     */
    while(!feof(fp))
    {
        char ichr[64];              //0.caution overflow.
        char isrc[64];              //1.
        char ityp[64];              //2.
        unsigned long ista;         //3.
        unsigned long iend;         //4.
        char isco[64];              //5.
        char iori[64];              //6.
        int ifra=0;                 //7. 0,1,2: ista+ifra meaning the real start of a CDS 
        char ides[1024];            //8.isoform should be derived from this item.
        char cfra[64];
        
        
        fscanf(fp, "%s\n",   ichr);
        if(ichr[0] == '#') continue;
        fscanf(fp, "%s\n",   isrc);
        fscanf(fp, "%s\n",   ityp);
        fscanf(fp, "%ld\n", &ista); // caution: the begining position of a real sequence is 1.
        fscanf(fp, "%ld\n", &iend); // caution.
        fscanf(fp, "%s\n",   isco);
        fscanf(fp, "%s\n",   iori);
        fscanf(fp, "%s\n",   cfra);
        fscanf(fp, "%s\n",   ides); 
        
        if(cfra[0] == '.') ifra = 0;
        else
        {
            ifra = cfra[0] - '0';
        }
        
        // caution: name consistence in all files
        if(reg_chromosome.compare((string)ichr)==0)
        if(iend>reg_begin && ista<reg_end) // as long as overlap found // caution on border
        {
            ////cout << ichr << "\t" << isrc << "\t" << ityp << "\t" << ista << "\t" << iend << "\t"
            ////<< isco << "\t" << iori << "\t" << ifra << "\t" << ides << endl;
            /* find gene name */
            /*
            std:string temp_str(ides);
            int nameB = temp_str.find("=")+1;
            int nameE = temp_str.find(";");  // "gene"
            if(nameE == std::string::npos)
            {
                nameE = temp_str.find(".");  // other
            }
            if(nameE == std::string::npos)
            {
                printf("NOT standard GFF format. Exit.\n");
                exit(1);
            }
            nameE = nameE - 1;
            std::string gene_name = temp_str.substr(nameB, nameE-nameB+1);
            */

            std:string temp_str(ides);
            
            int nameB = temp_str.find("=")+1;
            int nameE = temp_str.find(";");  // "gene"
            int nameE2= temp_str.find(".");  // other
            if(nameE == std::string::npos && nameE2 == std::string::npos)
            {
                printf("NOT standard GFF format. Exited (in read_GFFinfo.cpp).\n");
                exit(1);
            }
            if(nameE2 != std::string::npos && nameE != std::string::npos)
            {
                if(nameE2 < nameE)
                    nameE = nameE2;          // other
            }
            else if (nameE == std::string::npos)
            {
                nameE = nameE2;
            }  
            nameE = nameE - 1;
            std::string gene_name = temp_str.substr(nameB, nameE-nameB+1); // like "AT1G80865"
            
            /* set info of each seq_type */
            if(strncmp(ityp, "gene", 4)==0 || strncmp(ityp, "transposable_element_gene", 25)==0)  // caution: uppercase/lowercase of letters
            {                
                QUARTET gene_locus;
                gene_locus.start  = ista;
                gene_locus.end    = iend;
                gene_locus.orien += iori;
                gene_locus.seq   += chr_seq.substr(ista-1, iend-ista+1);// ista is a counted from 1; 
                                                                        // in c, string starts at 0.
                gene_ann.insert(std::pair<std::string, QUARTET>(gene_name, gene_locus));
                RANGE rge;
                 rge.start = ista;
                 rge.end   = iend;
                if(seq_type.find(gene_name) != seq_type.end())
                {
                    multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                    ann_itr = seq_type.find(gene_name);
                    (*ann_itr).second.insert(std::pair<std::string, RANGE>("gene", rge));
                }
                else
                {
                    multimap<std::string, RANGE> tmap;
                    tmap.insert(std::pair<std::string, RANGE>("gene", rge));
                    seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                               (gene_name, tmap));
                }
                ssnum = 0;
            }
            //else if(strncmp(ityp, "transcript", 10)==0) // for alpina: TODO; // isoform should be used: difference isofroms can have different mRNAs.
            else if(strncmp(ityp, "mRNA", 4)       ==0 || 
                    strncmp(ityp, "transcript", 10)==0 || 
                    strncmp(ityp, "ncRNA", 5)      ==0 || 
                    strncmp(ityp, "snoRNA", 6)     ==0 || 
                    strncmp(ityp, "tRNA", 4)       ==0 ||
                    strncmp(ityp, "rRNA", 4)       ==0 ||
                    strncmp(ityp, "miRNA", 5)      ==0 ||
                    strncmp(ityp, "snRNA", 5)      ==0 ||
                    strncmp(ityp, "mRNA_TE_gene", 12) == 0) // caution: non-coding RNA should be removed from analysis?
            {             
                if(RNAtype.find((string)ityp)==std::string::npos) 
                {
                    RNAtype += "#"+(string)ityp;  
                    if(verbose)
                    cout << "Warning: " << ityp << " related gene is included in analysis." << endl; 
                }
                std::string isoform = "";
                if(temp_str.find(";") != std::string::npos) // caution: gffs may have diff-formats.
                {
                    /* check: TAIR10_GFF3_gene.gff:        isoform=1,2,...   */
                    unsigned long lenform = temp_str.find(";")-1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                else 
                {
                    /* in alpina GFF, it should be 'transcript' seq_type???  */
                    /* check: A_alpina.V3.annotation.gff3: isoform=t1,t2,... */
                    unsigned long lenform = temp_str.length() -1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                
                RANGE rge;
                 rge.start = ista;
                 rge.end   = iend;
                std::string gene_name_isoform = gene_name + "." + isoform;                    // new

              //if(seq_type.find(gene_name) != seq_type.end())
                if(seq_type.find(gene_name_isoform) != seq_type.end())                        // new
                {
                    multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                  //ann_itr = seq_type.find(gene_name);
                    ann_itr = seq_type.find(gene_name_isoform);                               // new
                    (*ann_itr).second.insert(std::pair<std::string, RANGE>("mRNA", rge));
                    //printf("insert 2: mRNA/transcript of gene %s.\n", gene_name_isoform.c_str());
                }
                else
                {
                    multimap<std::string, RANGE> tmap;
                    tmap.insert(std::pair<std::string, RANGE>("mRNA", rge));
                  //seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                  //                           (gene_name, tmap));
                    seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                       (gene_name_isoform, tmap));                            // new
                    //printf("insert 1: mRNA/transcript of gene %s.\n", gene_name_isoform.c_str());
                }
                // update range of gene if necessary
                // assert(gene_ann.find(gene_name) != gene_ann.end());
                if(gene_ann.find(gene_name) == gene_ann.end())
                {
                    cout << "   Warning: no gene name found for " << gene_name << endl;
                    cout << ichr << "\t" << isrc << "\t" << ityp << "\t" << ista << "\t" << iend << "\t"
                         << isco << "\t" << iori << "\t" << ifra << "\t" << ides << endl;
                }
                
                bool changed = false;
                if(gene_ann.find(gene_name) != gene_ann.end() && ista < gene_ann[gene_name].start)
                {
                    gene_ann[gene_name].start = ista;
                    changed = true;
                    //cout << "   Info: gene " << gene_name << " updated as from " << ista << endl;                    
                }
                if(gene_ann.find(gene_name) != gene_ann.end() && iend > gene_ann[gene_name].end)
                {
                    gene_ann[gene_name].end = iend;
                    changed = true;
                    //cout << "   Info: gene " << gene_name << " updated as to " << iend << endl;                    
                }             
                if(changed)
                {
                    gene_ann[gene_name].seq   = chr_seq.substr(gene_ann[gene_name].start - 1, 
                                                               gene_ann[gene_name].end - gene_ann[gene_name].start + 1);
                                                               // ista is a counted from 1;
                }
            }
            else if(strncmp(ityp, "CDS", 3)==0) // isoform should be used
            {
                std::string isoform = "";
                if(temp_str.find(",") != std::string::npos) // caution: gffs may have diff-formats.
                {
                    /* check: TAIR10_GFF3_gene.gff:        isoform=1,2,...   */
                    unsigned long lenform = temp_str.find(",")-1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                else
                {
                    /* check: A_alpina.V3.annotation.gff3: isoform=t1,t2,... */
                    unsigned long lenform = temp_str.length() -1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                
                std::string gene_name_isoform = gene_name + "." + isoform;                    // new
                //if(isoform.find("1")==isoform.length()-1)   // caution: how if any cases like t11?
                //{            
                    /* set info of coding region */
                    QUARTET cds;
                    cds.start  = ista;
                    cds.end    = iend;
                    cds.orien += iori;
                    cds.seq   += chr_seq.substr(ista-1, iend-ista+1);
                    cds.frame  = ifra;
                    /* collect the above info in a coding map of gene_name */
                  //if(coding_ann.find(gene_name) != coding_ann.end())
                    if(coding_ann.find(gene_name_isoform) != coding_ann.end())                // new          
                    {
                        map<std::string, map<unsigned long, QUARTET> >::iterator coding_itr;
                      //coding_itr = coding_ann.find(gene_name);
                        coding_itr = coding_ann.find(gene_name_isoform);                      
                        (*coding_itr).second.insert(std::pair<unsigned long, QUARTET>(ista, cds));
                    }
                    else
                    {
                        map<unsigned long, QUARTET>  nullCMap;
                        nullCMap.insert(std::pair<unsigned long, QUARTET>(ista, cds));
                      //coding_ann.insert(std::pair<std::string, map<unsigned long, QUARTET> > 
                      //                             (gene_name, nullCMap));
                        coding_ann.insert(std::pair<std::string, map<unsigned long, QUARTET> > 
                                                     (gene_name_isoform, nullCMap));          // new                                            
                    }
                    /* collect seq type info in another map of gene_name */                    
                    RANGE rge;
                     rge.start = ista;
                     rge.end   = iend;
                    //if(seq_type.find(gene_name) != seq_type.end())
                    if(seq_type.find(gene_name_isoform) != seq_type.end())                    // new                    
                    {
                        multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                      //ann_itr = seq_type.find(gene_name);
                        ann_itr = seq_type.find(gene_name_isoform);                           // new
                        (*ann_itr).second.insert(std::pair<std::string, RANGE>("CDS", rge));
                    }
                    else
                    {
                        multimap<std::string, RANGE> nullSTMap;
                        nullSTMap.insert(std::pair<std::string, RANGE>("CDS", rge));
                      //seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                      //                           (gene_name, nullSTMap));
                        seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                                   (gene_name_isoform, nullSTMap));           // new              
                    }
                //}
            }
            else if(strncmp(ityp, "exon", 4)==0) // isoform should be used
            {                
                std::string isoform = "";
                if(temp_str.find(",") != std::string::npos) // caution: gffs may have diff-formats.
                {
                    /* check: TAIR10_GFF3_gene.gff:        isoform=1,2,...   */
                    unsigned long lenform = temp_str.find(",")-1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                else
                {
                    /* check: A_alpina.V3.annotation.gff3: isoform=t1,t2,... */
                    unsigned long lenform = temp_str.length() -1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                 
                // TODO: run on TAIR10_GFF2_genes.gff related data to see if condition below works. 
                // - 2013-04-03 20:56
                std::string gene_name_isoform = gene_name + "." + isoform;                    // new ::start from this line to modify
                // printf("exon of gene %s\n", gene_name_isoform.c_str());
                multimap<std::string, multimap<std::string, RANGE> >::iterator STitr;
                multimap<std::string, RANGE>::iterator MRitr;
                //STitr     = seq_type.find(gene_name);
                STitr     = seq_type.find(gene_name_isoform);                                 // new
                
                if(STitr == seq_type.end())
                {
                    printf("Error: gff is not in correct format; %s of seq type missing for checking mRNA. \n", gene_name_isoform.c_str());
                    printf("       Order of features should be: gene followed by mRNA by exon by CDS in GFF. Sorry you may need to reformat it.\n");
                    exit(1);
                }
                
                MRitr     = (*STitr).second.find("mRNA");   // question: why there is no "mRNA" info
                if(STitr != seq_type.end())                 // in file A_alpina_V3.annotation.gff3
                if(MRitr != (*STitr).second.end())
                //if(isoform.find("1")==isoform.length()-1) // caution: how if any cases like t11? - all isoforms should be used
                {
                    /* find out the start and end of a possible splice site */
                    RANGE erge;
                    unsigned long mstart = (*MRitr).second.start;
                    unsigned long mend   = (*MRitr).second.end;
                    if(ista != mstart)
                    {
                        ssnum ++;
                        erge.start = ista-2; // "A"
                        erge.end   = ista-1; // "G"
                        
                        /* insert this splice site info into map seq_type */
                        char ssc[512];
                        sprintf(ssc, "splice_site_change_%ld\0", ssnum);
                      //if(seq_type.find(gene_name) != seq_type.end())
                        if(seq_type.find(gene_name_isoform) != seq_type.end())                // new                  
                        {
                            multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                          //ann_itr = seq_type.find(gene_name);
                            ann_itr = seq_type.find(gene_name_isoform);                       // new
                          //(*ann_itr).second.insert(std::pair<std::string, RANGE>
                          //                         ((std::string)ssc, erge));
                          (*ann_itr).second.insert(std::pair<std::string, RANGE>
                                                     ((std::string)ssc, erge));               // new
                        }
                        else
                        {
                           multimap<std::string, RANGE> nullSTMap;
                           nullSTMap.insert(std::pair<std::string, RANGE>((std::string)ssc, erge));
                         //seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                         //                        (gene_name, nullSTMap));
                           seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                                   (gene_name_isoform, nullSTMap));           // new                                            
                        }
                    }
                    if(iend != mend)
                    {
                        ssnum ++;
                        erge.start = iend+1; // "G"
                        erge.end   = iend+2; // "T"
                        
                        /* insert this splice site info into map seq_type */
                        char ssc[512];
                        sprintf(ssc, "splice_site_change_%ld\0", ssnum);
                      //if(seq_type.find(gene_name) != seq_type.end())
                        if(seq_type.find(gene_name_isoform) != seq_type.end())                // new
                        {
                            multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                          //ann_itr = seq_type.find(gene_name);
                            ann_itr = seq_type.find(gene_name_isoform);                       // new                     
                            (*ann_itr).second.insert(std::pair<std::string, RANGE>
                                                     ((std::string)ssc, erge));
                        }
                        else
                        {
                           multimap<std::string, RANGE> nullSTMap;
                           nullSTMap.insert(std::pair<std::string, RANGE>((std::string)ssc, erge));
                         //seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                         //                        (gene_name, nullSTMap));
                           seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                                   (gene_name_isoform, nullSTMap));           // new                                      
                        }
                    }
                }
            }
            else if(strncmp(ityp, "five_prime_UTR", 14)==0||strncmp(ityp, "three_prime_UTR", 15)==0)
            {
                std::string isoform = "";
                if(temp_str.find(",") != std::string::npos) // caution: gffs may have diff-formats.
                {
                    /* check: TAIR10_GFF3_gene.gff:        isoform=1,2,...   */
                    unsigned long lenform = temp_str.find(",")-1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                else
                {
                    /* check: A_alpina.V3.annotation.gff3: isoform=t1,t2,... */
                    unsigned long lenform = temp_str.length() -1 - temp_str.find(".")-1 + 1;
                    isoform += temp_str.substr(temp_str.find(".")+1, lenform);
                    if((*isof).length() == 0) (*isof)     = isoform;
                }
                
                std::string gene_name_isoform = gene_name + "." + isoform;                    // new ::start from this line to modify
                RANGE rge;
                 rge.start = ista;
                 rge.end   = iend;
              //if(isoform.find("1")==isoform.length()-1)   // caution: how if any cases like t11?
              //if(seq_type.find(gene_name) != seq_type.end())
                if(seq_type.find(gene_name_isoform) != seq_type.end())
                {
                    multimap<std::string, multimap<std::string, RANGE> >::iterator ann_itr;
                  //ann_itr = seq_type.find(gene_name);
                    ann_itr = seq_type.find(gene_name_isoform);                               // new      
                    (*ann_itr).second.insert(std::pair<std::string, RANGE>(ityp, rge));
                }
                else
                {
                    multimap<std::string, RANGE> nullSTMap;
                    nullSTMap.insert(std::pair<std::string, RANGE>(ityp, rge));
                  //seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                  //                               (gene_name, nullSTMap));
                    seq_type.insert(std::pair<std::string, multimap<std::string, RANGE> >
                                                   (gene_name_isoform, nullSTMap));           // new                                     
                }
            }
            else ; // protein, tRNA:                       check file TAIR10_GFF3_genes.gff; 
                   // transcript, start_codon, stop_codon: check file A_alpina_V3.annotation.gff3
        }
    }
    fclose(fp);
    
    /* check gene_ann according to seq_type info: are there more than isoforms? */
    map<std::string, QUARTET> gene_ann_tmp;
    map<std::string, QUARTET>::iterator gene_itr;
    map<std::string, QUARTET>::iterator gene_itr_end;
    gene_ann_tmp.clear();
    gene_itr     = gene_ann.begin();
    gene_itr_end = gene_ann.end();
    while(gene_itr != gene_itr_end)
    {
         multimap<std::string, multimap<std::string, RANGE> >::iterator st_itr;        //seq_type
         multimap<std::string, multimap<std::string, RANGE> >::iterator st_itr_end;
         multimap<std::string, multimap<std::string, RANGE> >::iterator st_itr_isf;    //isoform itr
         st_itr     = seq_type.find((*gene_itr).first);
         st_itr_end = seq_type.end();
         if(st_itr != st_itr_end)
         {
             multimap<std::string, multimap<std::string, RANGE> >::iterator st_itr_tmp;// caution: not checked here
             st_itr_tmp =   st_itr;
             st_itr_isf = ++st_itr;
             while(st_itr_isf != st_itr_end)
             {
                  unsigned long isofpos = (*st_itr_isf).first.find(".")+1;
                  unsigned long lengeid = isofpos - 1 - 1 - 0 + 1;
                  std::string   gene0id = (*st_itr_isf).first.substr(0, lengeid);
                  if((*gene_itr).first == gene0id)
                  {
                      gene_ann_tmp.insert(std::pair<std::string, QUARTET>((*st_itr_isf).first, (*gene_itr).second));
                      //gene_ann.insert(std::pair<std::string, QUARTET>((*st_itr_isf).first, (*gene_itr).second));
                  }
                  st_itr_isf ++;
             }
             //seq_type.erase(st_itr_tmp); // seq_type contains "CDS" "mRNA" only (no "gene")
         } 
         else
         {
             printf("%s not found.\n", (*gene_itr).first.c_str());
         }
         gene_itr ++;
    }
    gene_ann.clear();
    gene_ann.insert(gene_ann_tmp.begin(), gene_ann_tmp.end()); // TODO: is there a more efficient way?
    
    return true;
}
