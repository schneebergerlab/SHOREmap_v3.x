#include    <string>
#include    <vector>
#include       <map>
#include   <fstream>
#include   <sstream>
#include  <iostream>
#include <algorithm>
#include   <ctype.h>

#include "split_string.h"
#include "globals.h"

using namespace std;

// find number of indels and remove indel flags from the original base string
int find_indel_length(string str, map<int, int>* indel, string* str2);

char determine_mut_base(int nA, int nC, int nG, int nT, int nN, int nD);

int convert_shore_pileup(int runid, char* shore_pileup_file)
{
    // open pileup file
    std::ifstream fp(shore_pileup_file);
    if(!fp.is_open())
    {
        cout << "Cannot open file " << shore_pileup_file << " to read consensus. Exited." << endl;
        return 1;
    }
    
    // initialize a file where converted consensus info will be recorded.
    fstream outfp;
    std::stringstream ss("");
    ss << out_folder << runid << "_converted_consensus_summary.txt";
    outfp.open((ss.str()).c_str(), ios::out);
    if(!outfp.is_open())
    {
        cout << "Cannot open file " << ss.str() << " to write converted consensus info. Exited.\n";
        return 1;
    }
    
    while(fp.good())
    {
        std::string line("");
        getline(fp, line);
        
        if(line.size()==0 || line.find("#")!=std::string::npos)
        {
            continue;
        }
        
        vector<string> lineinfo = split_string(line, '\t');
        
        // caution: only major allele
        
        outfp << lineinfo[0] << "\t" << lineinfo[1] << "\t"; // 0.chr, 1.pos
        
        // find the base string after removing indels
        map<int, int> indel;
        string cleanBaseStr;
        find_indel_length(lineinfo[9], &indel, &cleanBaseStr);
        
        //outfp << lineinfo[9] << "\t " << cleanBaseStr << "\t";
        
        size_t nf=0, nr=0, raw_cov=0;
        
        nf = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), '.');
        nr = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), ',');
        
        size_t na=0, nc=0, ng=0, nt=0, nn=0;
        size_t nA=0, nC=0, nG=0, nT=0, nN=0, nD=0, nI=0;
        na = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'a');
        nA = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'A');
        
        nc = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'c');
        nC = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'C');
        
        ng = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'g');
        nG = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'G');
        
        nt = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 't');
        nT = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'T');
        
        nn = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'n');
        nN = std::count(cleanBaseStr.begin(), cleanBaseStr.end(), 'N');
        
        char refbase = lineinfo[2][0];
        if(refbase == 'A')      {na = nr; nA = nf;} 
        else if(refbase == 'C') {nc = nr; nC = nf;} 
        else if(refbase == 'G') {ng = nr; nG = nf;} 
        else if(refbase == 'T') {nt = nr; nT = nf;} 
        else ;//cout << "Warning: non-determined reference base happened at chr " << lineinfo[0] << ":" << lineinfo[1] << endl;
        
        map<int, int>::iterator itr;
        map<int, int>::iterator itr_end;
        itr     = indel.begin();
        itr_end = indel.end();
        while(itr != itr_end)
        {
            if((*itr).first < 0)
            {
                nD ++;
            }
            else
            {
                nI ++;
            }
            itr ++;
        }
        
        // determine base
        char base = determine_mut_base(na+nA, nc+nC, nG+ng, nT+nt, nn+nN, nD);
        
        // determine raw cov
        raw_cov = na+nA + nc+nC + nG+ng + nT+nt + nn+nN + nD;
        
        // 2. mut-base 3.raw-cov 4.A, 5.C, 6.G, 7.T, 8.Del, 9.N, 
        outfp <<  base  << "\t" << raw_cov << "\t"
              << na+nA  << "\t" << nc+nC << "\t" << ng+nG   << "\t" << nt+nT << "\t"
              <<  nD    << "\t" << nn+nN << "\t";
        // 10. avg_hit 11. avg_mm
        outfp << lineinfo[7] << "\t" << lineinfo[8] << endl;  
    }
    
    fp.close();
    outfp.close();
    return 0;
}

int find_indel_length(string str, map<int, int>* indel, string* str2)
{
    string str_indel_rmd("");
    for (size_t pos=0; pos<str.length(); pos++)
    {
        if(str[pos]=='+' || str[pos]=='-')
        {
            int tmpkey = pos+1; // caution: postion has become 1-based system
            char id = str[pos];
            pos ++;
            stringstream ss;
            while(isdigit(str[pos]))
            {
                ss << str[pos];
                pos++;
            }
            if(id=='-') tmpkey = tmpkey*(-1);
            int len = atoi(ss.str().c_str());
            (*indel).insert(std::pair<int, int>(tmpkey, len));            
            pos += len-1;
        }
        else
        {
            str_indel_rmd += str[pos];
        }
    }    
    *str2 = str_indel_rmd;
    return 0;
}

char determine_mut_base(int nA, int nC, int nG, int nT, int nN, int nD)
{
    char base = 'N';
    
    int  numlist[]  = {nA,  nC,  nG,  nT,  nN,  nD};
    char baselist[] = {'A', 'C', 'G', 'T', 'N', '-','R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'X'};
    
    int  nummax    = 0;
    int  mycount   = 0;
    char charmax[6];
    
    for(int i=0; i < 6; i++)
    {
        if(numlist[i] > nummax)
        {
            mycount   = 1;
            nummax    = numlist[i];
            charmax[mycount-1] = baselist[i];
        }
        else
        if(numlist[i] == nummax)
        {
            mycount   ++;
            charmax[mycount-1] = baselist[i];
        }
        else ;
    }
    if(mycount == 1)
    {
        base    = charmax[0];
    }
    else if(mycount == 2)
    {
        if(charmax[0]=='A' && charmax[1]=='G') base = baselist[6];  else
        if(charmax[0]=='C' && charmax[1]=='T') base = baselist[7];  else
        if(charmax[0]=='G' && charmax[1]=='T') base = baselist[8];  else
        if(charmax[0]=='A' && charmax[1]=='C') base = baselist[9];  else
        if(charmax[0]=='C' && charmax[1]=='G') base = baselist[10]; else
        if(charmax[0]=='A' && charmax[1]=='T') base = baselist[11]; else 
        base = 'X';
    }
    else if(mycount == 3)
    {
        if(charmax[0]=='C' && charmax[1]=='G' && charmax[2]=='T') base = baselist[12]; else
        if(charmax[0]=='A' && charmax[1]=='G' && charmax[2]=='T') base = baselist[13]; else
        if(charmax[0]=='A' && charmax[1]=='C' && charmax[2]=='T') base = baselist[14]; else
        if(charmax[0]=='A' && charmax[1]=='C' && charmax[2]=='G') base = baselist[15]; else
        base = 'X';
    }
    else if(mycount == 4)
    {
        base = 'N';
    }
    else
    {
        base = 'X';
    }
    
    return base;
}

// #chr	pos	refbase	GCcont	seqC	expcov	cov	avghits	avgmm	bases	qualities	mapT
// 1	25	T	-1	7	36.445	8	1	0	........	LLLLLLLL	11111113

/* IUPAC code for DNA
A	A
C	C
G	G
T	T
R	A or G
Y	C or T
K	G or T
M	A or C
S	C or G
W	A or T
B	not A (i.e. C, G, T)
D	not C (i.e. A, G, T)
H	not G (i.e. A, C, T)
V	not T (i.e. A, C, G)
N	A C G T
X	masked
-	gap of indeterminate length
*/











