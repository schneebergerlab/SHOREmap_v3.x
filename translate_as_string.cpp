/* this function translates a DNA sequence into protein sequence with standard coding table       */

#include <string>
#include <algorithm>
#include <map>
#include <stdio.h>

using namespace std;

void initialize_coding_table(map<std::string, std::string>* code);

void translate_as_string(std::string mydna, std::string* myprotein)
{
    /* initialize */
    map<std::string, std::string> transl_table;
    initialize_coding_table(&transl_table);
    
    (*myprotein).clear();
    for(unsigned long i = 0; i < mydna.size(); i += 3)
    {
        std::string codon = mydna.substr(i, 3);
        std::transform(codon.begin(), codon.end(),codon.begin(), ::toupper);
        
        if(transl_table.find(codon) != transl_table.end())
        {
            (*myprotein) += transl_table[codon];
        }
        else if (codon.length() == 3) // avoid null/false codon
        {
            (*myprotein) += "X"; // caution: unknown letters like 'N' appear in coding sequences
        }
        else ;
    }
}

/* this function initializes a table for translating a dna to protein                             */
/* TODO: check: are any other non-standard codes useful?                                          */
void initialize_coding_table(map<std::string, std::string>* code)
{
    /* standard code: 
       http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1 
    */
    std::string AAs    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    std::string Starts = "---M---------------M---------------M----------------------------";
    std::string Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    std::string Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    std::string Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
  
    for(int i = 0; i < AAs.size(); i ++)
    {
        std::string trip("");
        trip = Base1.substr(i,1) + Base2.substr(i,1) + Base3.substr(i,1);
        (*code).insert(std::pair<std::string, std::string>(trip, (std::string)AAs.substr(i,1)));
    }
}
