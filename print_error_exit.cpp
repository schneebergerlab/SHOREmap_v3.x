#include <stdio.h>
#include <stdlib.h>

////// if there is an error on an option, exit the program //////
void print_error_exit(char* opt, bool fscanf_error)
{
    if(fscanf_error)
    {
        printf("ERROR on fscanf in %s. Exited.\n", opt);
    }
    else // cmd options
    {
        printf("ERROR: %s requires (proper) argument. Exited.\n", opt);
    }
    exit(1);
}
