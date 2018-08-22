debug_noisy:
	g++ -g -o SHOREmap SHOREmap.cpp is_number.cpp init_backcross.cpp filter_fg.cpp ShoreMap_annotate.cpp plot_chr_snp_annotations.cpp ShoreMap_convert.cpp convert_shore_pileup.cpp globals.cpp compare2strings.cpp read_chromosomes.cpp read_referror.cpp read_marker.cpp read_allele_count.cpp write_log.cpp write_zoom_region.cpp print_allele_count.cpp read_allele_count2.cpp ShoreMap_extract.cpp precheck_opt.cpp print_error_exit.cpp read_extract_cons.cpp read_extract_bg_qref.cpp ShoreMap_outcross.cpp cmd_initialize.cpp filter_with_marker_parent.cpp plot_chr_winboost.cpp print_plot_info.cpp init_SNP.cpp get_SNPlist.cpp read_chromosome_seq.cpp read_GFFinfo.cpp translate_as_string.cpp get_gene_snps.cpp get_protein_changes.cpp get_INDELlist.cpp init_INDEL.cpp check_ref_err.cpp split_string.cpp print_filtered_marker.cpp read_peaks.cpp ShoreMap_create2.cpp find_marker_pos_2parents_v3.cpp kmeans.cpp Calc_Wilson_score_interval.cpp find_z_value.cpp  ShoreMap_idFilter.cpp -L/usr/lib/ -lXt -L./dislin -ldislin_d -lm
	
	#-B -g -DVERBOSE_DEBUG -DDEBUG

release:
	g++ -g -o SHOREmap SHOREmap.cpp is_number.cpp init_backcross.cpp filter_fg.cpp ShoreMap_annotate.cpp plot_chr_snp_annotations.cpp ShoreMap_convert.cpp convert_shore_pileup.cpp globals.cpp compare2strings.cpp read_chromosomes.cpp read_referror.cpp read_marker.cpp read_allele_count.cpp write_log.cpp write_zoom_region.cpp print_allele_count.cpp read_allele_count2.cpp ShoreMap_extract.cpp precheck_opt.cpp print_error_exit.cpp read_extract_cons.cpp read_extract_bg_qref.cpp ShoreMap_outcross.cpp cmd_initialize.cpp filter_with_marker_parent.cpp plot_chr_winboost.cpp init_SNP.cpp get_SNPlist.cpp read_chromosome_seq.cpp read_GFFinfo.cpp translate_as_string.cpp get_gene_snps.cpp get_protein_changes.cpp get_INDELlist.cpp init_INDEL.cpp check_ref_err.cpp split_string.cpp print_filtered_marker.cpp read_peaks.cpp ShoreMap_create.cpp find_marker_pos_2parents.cpp  Calc_Wilson_score_interval.cpp find_z_value.cpp  -lXt -L./dislin10.3 -ldislin_d -lm
	
	
## static 2014-07-30 -- check src/SHOREmap_static_test

# with removed files:
#debug_noisy:
#	g++ -B -g -DVERBOSE_DEBUG -DDEBUG -o SHOREmap SHOREmap.cpp is_number.cpp init_outcross.cpp init_backcross.cpp filter_fg.cpp ShoreMap_annotate.cpp ShoreMap_convert.cpp convert_shore_pileup.cpp globals.cpp read_chromosomes.cpp read_referror.cpp read_marker.cpp read_allele_count.cpp write_log.cpp write_zoom_region.cpp print_allele_count.cpp SHOREmap_plot.cpp ShoreMap_confint.cpp filterSampling.cpp windowFinding.cpp calc_region_freqs.cpp find_max_Mod.cpp find_cur_window.cpp calc_pbinom.cpp calc_windows_avg.cpp identify_peaks.cpp multll.cpp extend_interval.cpp maxConf.cpp restrictedModel.cpp index_dataset.cpp read_allele_count2.cpp ShoreMap_extract.cpp precheck_opt.cpp print_error_exit.cpp read_extract_cons.cpp read_extract_bg_qref.cpp ShoreMap_visualize.cpp cmd_initialize.cpp plot_chr_AF.cpp plot_chr_winboost.cpp init_SNP.cpp get_SNPlist.cpp read_chromosome_seq.cpp read_GFFinfo.cpp translate_as_string.cpp get_gene_snps.cpp get_protein_changes.cpp get_INDELlist.cpp init_INDEL.cpp check_ref_err.cpp split_string.cpp  Calc_Wilson_score_interval.cpp find_z_value.cpp -lXt -ldislin -lm

# removed 16 files
# ShoreMap_confint.cpp
# identify_peaks.cpp 
# multll.cpp 
# extend_interval.cpp 
# maxConf.cpp 
# restrictedModel.cpp 
# filterSampling.cpp
# find_max_Mod.cpp 
# find_cur_window.cpp 
# calc_pbinom.cpp
# SHOREmap_plot.cpp
# windowFinding.cpp
# calc_windows_avg.cpp
# calc_region_freqs.cpp
# index_dataset.cpp
# plot_chr_AF.cpp (this is replace with plot_chr_winboost.cpp)

#before 2013-09-17 10:10am

#debug_noisy:
#	g++ -g -o SHOREmap SHOREmap.cpp is_number.cpp init_backcross.cpp filter_fg.cpp ShoreMap_annotate.cpp ShoreMap_convert.cpp convert_shore_pileup.cpp globals.cpp compare2strings.cpp read_chromosomes.cpp read_referror.cpp read_marker.cpp read_allele_count.cpp write_log.cpp write_zoom_region.cpp print_allele_count.cpp read_allele_count2.cpp ShoreMap_extract.cpp precheck_opt.cpp print_error_exit.cpp read_extract_cons.cpp read_extract_bg_qref.cpp ShoreMap_outcross.cpp cmd_initialize.cpp filter_with_marker_parent.cpp plot_chr_winboost.cpp init_SNP.cpp get_SNPlist.cpp read_chromosome_seq.cpp read_GFFinfo.cpp translate_as_string.cpp get_gene_snps.cpp get_protein_changes.cpp get_INDELlist.cpp init_INDEL.cpp check_ref_err.cpp split_string.cpp print_filtered_marker.cpp read_peaks.cpp ShoreMap_create.cpp find_marker_pos_2parents.cpp  Calc_Wilson_score_interval.cpp find_z_value.cpp -L/usr/lib/ -lXt -L./dislin10.3 -ldislin_d -lm

#2014-03-28 19:34
# find_marker_pos_2parents_v2.cpp is replaced by find_marker_pos_2parents_v3.cpp
