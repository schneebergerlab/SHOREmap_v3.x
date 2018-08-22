// pls check the header file: read_chromosome_seq.h for more info about this function. //
#include   <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include    <string>
#include   <fstream>
#include   <iostream>
#include  <assert.h>

using namespace std;

bool read_chromosome_seq(char* fgenome, char* chr_id, unsigned long max_chr_len, std::string* chr_seq)
{
    /* max_chr_len is not used - TO REMOVE - 2013-05-17 15:02                                     */
    
    std::ifstream fp (fgenome);
    if(!fp.is_open())
    {
        printf("   Error: reference genome file \'%s\' does NOT exist. Exited.\n", fgenome);
        exit(1);
    }
    
    bool ifound = false;
    while(fp.good())
    {
        std::string line("");
        getline(fp, line);
        if(line.size()==0) continue;
        
        size_t pos = line.find(" ");
        string seqname = line.substr(1);
        
        if(line[0]=='>' && pos!=std::string::npos)
        {
            printf("   Warning: you may have redundant info after \'>\' as a sequence name; need to remove. ");
            assert(pos>=2);
            seqname = line.substr(1, pos-1);
        }
        
        if(seqname.compare((std::string)chr_id) == 0) // after '>' it is chromosome id
        {
            //cout << "   Info: seq id " << seqname << " matched. Getting sequence..." << endl;
            line.clear();
            getline(fp, line);
            while(line.find(">") == std::string::npos && fp.good())
            {
                (*chr_seq) += line;
                line.clear();
                getline(fp, line);
            }
            ifound = true;
            break;
        }
    } 
    fp.close();
    if(!ifound) 
    {
        printf("   Error: Could not find genome sequence under given id %s. Exited!\n.", chr_id);
        printf("   Hint: please keep seq-ids (without whitespaces) in these four files of markers, fasta, chrSizes, and GFF exactly the same.\n");
        exit(1);
    }
    else        return  true;
}

/* caution: chr_id should be the same as in the genome file; otherwise error can happen           */
/* TODO: how to make it become a standard format???                                               */
