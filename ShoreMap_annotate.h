/* this function annotates given SNPs(, Indels) according to gene annotation file */
/*
Usage: SHOREmap annotate 

Mandatory:
--snp      STRING     Point mutation file
--chrom    STRING     Chromosome of target region
--start    INT        Start of target region
--end      INT        End of target region

Functional annotation:
--genome   STRING     Reference sequence file 
--gff      STRING     Gene annotation in GFF format

Optional:
--referr   STRING     Reference sequence errors
--del      STRING     Deletion file  // TODO 2013-03-26
--ins      STRING     Insertion file // TODO 2013-03-26
*/

void ShoreMap_annotate(int argc, char* argv[]);


/*
    std::string   ecotype = "";
    std::string   chromosome = "";
    unsigned long position = 0;
    std::string   stype = "intergeneric";
    std::string   gene_id = "NA";
    unsigned long gene_pos=0;
    unsigned long cds_pos=0;
    unsigned long codon_pos = 0;
    unsigned long ns_change = 0;
    unsigned long new_stop = 0;
    unsigned long lost_stop = 0;
    unsigned long splicechange = 0;
    std::string   ref_base="";
    std::string   new_base="";
    std::string   ref_aa ="";
    std::string   new_aa="";
    //map<std::string, std::string>domain_change;
    unsigned long support = 0;
    double        concordance = 0;
    double        quality = 0;
    unsigned long peak_distance = 999999999;
    double        marker_ratio = 0;
    std::string   source = "IlluminaGA2";

*/
