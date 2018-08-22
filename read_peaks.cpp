/* this function reads in peaks recorded in plot_chr_winboost.cpp step */

#include      <map>
#include   <string>
#include   <vector>
#include  <fstream>
#include  <stdio.h>
#include <stdlib.h>

#include "split_string.h"

using namespace std;

bool read_peaks(char* peak_file, multimap<std::string, unsigned long>* mpeaks)
{
    std::ifstream fpeaks (peak_file);
    if(!fpeaks.is_open())
    {
        printf("Cannot open file %s to read peaks (in read_peaks.cpp). ", peak_file);
        printf("No peaks recorded (thus 1 used as default). \n");
        return false;
    }
    std::string line;
    while(fpeaks.good())
    {
        getline(fpeaks, line);
        if (line.length() == 0) {continue;} // null line
        vector<string> infoline = split_string(line, '\t');
        if(infoline.size() != 2) continue;

        (*mpeaks).insert(std::pair<std::string, unsigned long>(infoline[0], atol(infoline[1].c_str())));
    }
    fpeaks.close();
    return true;
}
