/* read gene info from the annotation file                                                        */
/* GFF format (nine columes for each item):                                                       */
/*     0.chr 1.source 2.seq_type 3.start 4.end 5.score 6.orientation 7.frame 8.description        */

bool read_GFFinfo(char* file_gff, std::string* isof);
