#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include "print_error_exit.h"
#include "is_number.h"

#include "globals.h"
//////  pre-check if an option is valid; if yes, put it in map //////
int precheck_opt(int argc, char* argv[], int opt_i, char* opt_name, bool arg, bool is_num)
{
    if(arg)
    {
        opt_i ++;
        if (opt_i >= argc) print_error_exit(opt_name, false);
        if(is_num)
        {
            if(!is_number(argv[opt_i])) {
                printf("%s should be a number. Exited.\n", opt_name); 
                exit(1);
            }
        }
        // record <option, arg> in map<string, string>//
        string s_key(opt_name); 
        string s_val(argv[opt_i]);
        CMD.insert(std::pair<string,string>(s_key, s_val));
    }
    else
    {
        // record <option, "true"> in map<string, string>//
        string s_key(opt_name); 
        CMD.insert(std::pair<string,string>(s_key, "true"));
    }

    return opt_i;
}
