#include <map>

bool find_marker_pos_2parents_v3(char* paSNPfile, 
                                 char* pbREFfile, 
                                 char* pbSNPfile,
                                 char* paREFfile,
                                 char* F2SNPfile,
                                 unsigned long min_cov,
                                 unsigned long max_cov,
                                 double        min_ccd, 
                                 unsigned long max_N, 
                                 unsigned long max_ID, 
                                 unsigned long min_refq,
                                 std::string*  F2marker_ParentAsumBminusAB);
