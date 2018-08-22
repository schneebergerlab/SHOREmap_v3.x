/* this function records filtered markers such as before visualization if required from cmd       
   inputs:
      consen_info: e.g., CHR2POS2_ale1_ale2_err_COUNT or its subset from read_allele_counts2(...)
      ref_allele : e.g., ALLELE2                      or its subset from read_marker(...)
      mut_allele : e.g., ALLELE1                      or its subset ...
      mut_quality: e.g., QUALITY1                     or its subset ...
   Date 2013-05-28                                                                                */

bool print_filtered_marker(map<std::string, map<unsigned long, TRIPLE> > consen_info,
                           map<std::string, std::string> ref_allele,
                           map<std::string, std::string> mut_allele,
                           map<std::string, std::string> mut_quality);
