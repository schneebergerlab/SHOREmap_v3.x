/* input:  file contains ids and sequences of chromosomes */
/* output: <chr_id, string> in chr_seq                    */
/* date:   2013-04-02                                     */
/* file format:
               >chr_id1
               sequence1
               >chr_id2
               sequence2
               .
               .
               .
*/
#include <string>
bool read_chromosome_seq(char* fgenome, char* chr_id, unsigned long max_chr_len, std::string* chr_seq);
