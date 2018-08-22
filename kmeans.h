/***** The kmeans used in SHOREmap is modified by  SHQ from  the  origial  implementation
** kmeans.c
** - a  simple  k-means  clustering  routine  of  1-dimension  data  of  k-mer  frequency
** - returns the (minimum, maximum) paired-values of clusters of the data points in a map
**   
** Parameters
** - array of data points (double **data)
** - number  of  data  points   (int   n)
** - dimension        (int             m)
** - desired number of clusters  (int  k)
** - error     tolerance    (double    t)
**   - used as the  stopping  criterion, i.e. when  the  sum  of
**     squared euclidean distance  (standard  error for k-means)
**     of an iteration is within the  tolerable  range from that
**     of the previous iteration, the  clusters  are  considered
**     "stable", and the function returns
**   - a suggested value would be 10
** - output address for the final centroids (double **centroids)
**   - user must make sure the memory is properly allocated,  or
**     pass the null pointer if not interested in the  centroids
** Notes
** - this  function  is  provided  as  is  with   no   warranty.
** - the  author  is  not  responsible  for  any  damage  caused
**   either directly  or  indirectly  by  using  this  function.
** - anybody is free to  do  whatever  he/she  wants  with  this
**   function as long  as  this  header  section  is  preserved.
** Created on 2005-04-12 by  -  Roger  Zhang  (rogerz@cs.dal.ca)
** 
** Last compiled under g++ (Ubuntu/Linaro 4.7.2-2ubuntu1)  4.7.2
*/

int k_means(double **data, 
            unsigned long n, 
            unsigned long m, 
            unsigned long k, 
            double t, 
            double **centroids,
            unsigned long* cluster_counts,
            double *myScorSet);
