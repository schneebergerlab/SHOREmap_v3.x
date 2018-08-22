/* this function looks for a z value corresponding to a confidence percentile 

   date: 2015-06-15 15:04
*/

#include <math.h>

double find_z_value(double z_table[32][11], double confidence_percentile)
{
    double percentile_z = 1.96; // default percentile_z for confidence 0.95
    int    row     = 0;
    int    col     = 0;
    double minimum = 1;
    
    for(int i=0; i<31; i++)     // last row is not z percentile
    {
        for(int j=1; j<11; j++) // first column is z percentile
        {
            if(fabs(z_table[i][j]-confidence_percentile) < minimum)
            {
                row = i;
                col = j;
                minimum = fabs(z_table[i][j]-confidence_percentile);
            }
        }
    }
    percentile_z = z_table[row][0] + z_table[31][col];
    return percentile_z;
}
