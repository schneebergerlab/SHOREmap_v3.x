/* input  : file contains info about SNPs of a parent/background (e.g., quality_variant.txt)
   funtion: filter mutant-SNPs with SNPs contained in the given file                     
            case 1: this parent X          shows the same phenotype as the mutant;
            case 2: this parent X does not show  the same phenotype as the mutant.
            for case 1: remove SNPs in the mutant but not in X;
            for case 2: remove SNPs in the mutant as well as X.
            note: SNPs of parents can be controlled by using its mut-base qual, ccd, cov, etc.
   date   : 2013-May-24
   file format: project_name  chr_id  position ref_base allele_base other1(quality, coverage, etc)
                 
                E.G.:
                 
                proj_1        chr1    1234567  A        C	    ignored...
                proj_1        chr2    4567123  C        G           ignored...
                proj_1        chr3    6712345  G        T           ignored...
                proj_1        chr4    3456712  T        A           ignored...
                proj_1        chr5    2345671  N        C           ignored...      
                                          note:N
   note1: at least 5 columns while only last 4 columns from chr to allele_base are recorded.
   note2: SNPs of this parent should have high quality, coverage(, frequency) etc.
   note3: mutant-SNPs have been read into global CHR2POS2_ale1_ale2_err_COUNT beforehand.
*/

bool filter_with_marker_parent(char* fmarker_p, bool iphenotype);
