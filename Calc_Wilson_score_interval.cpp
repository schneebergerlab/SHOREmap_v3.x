/* this function calculates the confidence interval for the proportion of mutant base coverage 
   
   date: 2015-06-15 15:02
*/

#include <math.h>

int Calc_Wilson_score_interval(int n, int mut_cov, double percentile_z, double interval[2])
{
    // n           : sample size
    //               = number of reads of covering mutant and reference alleles at a marker position
    // mut_cov     : number of reads covering that marker position
    // percentile_z: 1-alpha/2 percentile of a standard normal distribution, alpha being error percentile
    //               percentile_z and alpha can be checked from a known table
    // estimated_p : is the proportion of successes in a Bernoulli trial process estimated from sample
    //               = ratio of mutant allele
    
    double square_z    = pow(percentile_z, 2);
    double factor_a    = 1/(1+square_z/n);
    double estimated_p = (double)mut_cov/n;
    
    double square_n    = (double)pow(n, 2);
    
    double factor_s    = sqrt(estimated_p*(1-estimated_p)/n + square_z/(4*square_n));
    
    interval[0]        = factor_a*(estimated_p + square_z/(2*n) - percentile_z*factor_s);
    interval[1]        = factor_a*(estimated_p + square_z/(2*n) + percentile_z*factor_s);
    
    return 0;
}
