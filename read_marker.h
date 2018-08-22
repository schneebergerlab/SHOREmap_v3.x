/* input : file contains info about markers (e.g., quality_variant.txt from SHORE)                */
/* output: <chr.pos, allele> in global variable ALLELE1 and ALLELE2                               */
/*       : number of markers given in the file, and number of markers to use                      */
/* date  : 2013-March-07                                                                          */
/* 
   file format: project_name  chr_id  position ref_base allele_base other1(quality, coverage, etc)
                 
                E.G.:
                 
                proj_1        chr1    1234567  A        C	     ignored...
                proj_1        chr2    4567123  C        G           ignored...
                proj_1        chr3    6712345  G        T           ignored...
                proj_1        chr4    3456712  T        A           ignored...
                proj_1        chr5    2345671  N        C           ignored...      
                                          note:N
   note: at least given 5 columns but only 4 columns from chr to allele_base are recorded.
*/
bool read_marker(char* fmarker, unsigned long* num_given, unsigned long* num_inuse);
