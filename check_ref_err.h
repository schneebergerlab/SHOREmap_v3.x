/* this function filters targeted markers of one mutant according to another mutant based on below*/
/*   With two sequenced mutant genomes, we use the targeted one as foreground, and the other as   */
/*   background to reduce the number of markers for analysis. Specifically, suppose we have a     */
/*   marker M in mutant A, we check if M also appears in mutant B. If it appears in B, we could   */
/*   remove M for future analysis as if it were resulted from an error residing in the reference  */
/*   genome. If it does not appear in mutant B, we further check quality_reference of mutant B.   */
/*   If the reference allele of marker M appears in quality_reference of mutant B but with low    */
/*   support/concordance, we might discard this marker M.                                         */
/*   Otherwise, we keep marker M for future analysis.                                             */
/* Date 2013-04-16 */

int check_ref_err(char*  fg_snp_file,
                  char*  fg_consencall_file, 
                  char*  bg_ref_base_file, 
                  unsigned long min_cov,
                  unsigned long max_cov,
                  double min_ccd, 
                  unsigned long max_N, 
                  unsigned long max_ID, 
                  unsigned long min_refq,
                  std::string* fg_snp_file_filtered,
                  map<std::string, MARK6> * fg_snp_map_filtered);
