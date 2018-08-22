/* this function plots boost value (of window-markers) along a chromosome                         */
/* chrID is the real id of a chromosome to plot, while plot_ith_chr indicates how many chromosomes*/
/* to plot if the current one included                                                            */
/* Date  : 2013-April-24                                                                          */

bool plot_chr_winboost(std::string chrID, 
                       unsigned long plot_ith_chr,
                       map<unsigned long, TRIPLE> internalData, 
                       map<unsigned long, TRIPLE> filtered);
