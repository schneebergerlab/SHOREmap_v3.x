/* for more info, please check the header file kmeans.h                                           */

#include <stdlib.h>
#include  <stdio.h>
#include <assert.h>
#include  <float.h>
#include   <math.h>
#include      <map>

#define FMAX 999999999.0

using namespace std;

int k_means(double **data, 
            unsigned long n, 
            unsigned long m, 
            unsigned long k, 
            double t, 
            double **centroids,
            unsigned long* cluster_counts,
            double *myScorSet)
{  
   printf("\tStarting clustering...\n");
   /* initialization                                                                              */
   unsigned long h, i, j;                                                  /* loop       counters */
   unsigned long *labels = (unsigned long*)malloc(n*sizeof(unsigned long));/*label for each point */
   if(!labels)
   {
       printf("Cannot allocate memory for labels with n=%ld\n", n);
       exit(1);
   }
   double old_error, error = FMAX;                           /* sum of squared euclidean distance */
   double **c1 = (double**)malloc(k*sizeof(double*));                      /* temp      centroids */
   if(!c1)
   {
       printf("Cannot allocate memory for c1 with k=%ld\n", k);
       exit(1);
   }
   assert(data && k > 0 && k <= n && m > 0 && t >= 0);                     /* for       debugging */
   
   for (h = 0, i = 0; i < k; h += n / k, i++) 
   {
      c1[i] = (double*)malloc(m*sizeof(double));
      assert(c1[i]);
      /* initialize k points as centroids                                                         */
      for (j = m; j-- > 0; centroids[i][j] = data[h][j]);
      ////for (j = 0; j< m; j++)
      ////{
      ////    printf("\t%.2f ", centroids[i][j]);
      ////}
      ////printf("\n");
   }
   /* loop of clustering                                                                          */
   unsigned long num_iteration = 0;
   do {
          old_error = error, error = 0;                                    /* update        error */
          /* clear old counts and temp centroids                                                  */
          for (i = 0; i < k; cluster_counts[i++] = 0) 
          {
              for (j = 0; j < m; c1[i][j++] = 0);
          }
          for (h = 0; h < n; h++) 
          {
              /* identify the closest cluster                                                     */
              double min_distance = DBL_MAX;
              for (i = 0; i < k; i++) 
              {
                 double distance = 0;
                 for (j = m; j-- > 0; distance += pow(data[h][j] - centroids[i][j], 2)) ;
                 if (distance < min_distance) 
                 {
                    labels[h] = i;
                    min_distance = distance;
                 }
              }
              /* update size and temp centroid of the destination cluster                         */
              for (j = m; j-- > 0; c1[labels[h]][j] += data[h][j]);
              cluster_counts[labels[h]]++;
              /* update standard error                                                            */
              error += min_distance;
           }

           for (i = 0; i < k; i++) 
           {   /* update all centroids                                                            */
               for (j = 0; j < m; j++) 
               {
                   centroids[i][j] = cluster_counts[i] ? c1[i][j] / cluster_counts[i] : c1[i][j];
               }
           }
           num_iteration ++;
           if(num_iteration > (unsigned long)FMAX)
           {
                printf("Max number of iterations reached.\n");
                break;
           } 
   } while (fabs(error - old_error) > t);
   
   /* find the maximum values of each attribute (cov, base-quality etc)  for each cluster so that we 
      can show it on legend... -TODO 2014-03-18 10:33                                             */
   
   
   /* rank clusters according element-summation of cluster centers                                */
   double* cluster_score = (double*)malloc((k+1)*sizeof(double));
   int*    cluster_id    = (int*)malloc((k+1)*sizeof(int));
   assert(cluster_score && cluster_id);
   //printf("\t cluster centers before sorting: \n");
   for(int ki=0; ki<k; ki ++)
   {
        double kscore  = 0.0;
        double mkscore = 1.0;
        //printf("\t cluster %d: ", ki);
        for(int cenj=0; cenj<m; cenj++)
        {
            //printf("%.2f ", centroids[ki][cenj]);
            kscore  += centroids[ki][cenj];
            mkscore *= centroids[ki][cenj];
        }
        cluster_score[ki] = kscore;
        cluster_id[ki]    = ki;
        //printf("old: cluster %d score %.2f count %d\n", 
        //             cluster_id[ki], cluster_score[ki], cluster_counts[ki]);
   }
   
   /* sort scores */
   double* tmp_center = (double*)malloc(m*sizeof(double));
   assert(tmp_center);
   for(int si=0; si<k; si++)
   {
       for(int sj=si+1; sj<k; sj++)
       {
           if(cluster_score[sj] < cluster_score[si]) // swap
           {
               // swap cluster_score
               double tmp_score  = cluster_score[si];
               cluster_score[si] = cluster_score[sj];
               cluster_score[sj] = tmp_score;
               // swap cluster_id
               int    tmp_cid    = cluster_id[si];
               cluster_id[si]    = cluster_id[sj];
               cluster_id[sj]    = tmp_cid;
               // swap counts
               int    tmp_cnt    = cluster_counts[si];
               cluster_counts[si]= cluster_counts[sj];
               cluster_counts[sj]= tmp_cnt;
               // swap centroids
               for(int tmpci=0; tmpci<m; tmpci++)
               {
                   double tmp_element   = centroids[si][tmpci];
                   centroids[si][tmpci] = centroids[sj][tmpci];
                   centroids[sj][tmpci] = tmp_element;
               } 
           }
       }
   }
   
   map<int, double> cluster_score_new;
   for(int si=0; si<k; si++)
   {
       cluster_score_new.insert(std::pair<int, double>(cluster_id[si], (double)si));
   }
   
   /*
   printf("\t cluster centers after sorting: \n");
   for(int ki=0; ki<k; ki ++)
   {
        printf("\t cluster %d: ", ki);
        double mkscore = 1.0;
        for(int cenj=0; cenj<m; cenj++)
        {
            printf("%.2f ", centroids[ki][cenj]);
        }
        printf(" swapped.double=%.2f new.int=%.0f count=%d\n",  
                 cluster_score[ki], cluster_score_new[cluster_id[ki]], cluster_counts[ki]);
   }
   */
   
   for(unsigned long di=0; di<n; di++)
   {
       int cid       = labels[di];
       map<int, double>::iterator csn_itr;
       csn_itr       = cluster_score_new.find(cid);
       myScorSet[di] = (*csn_itr).second; 
   }

   /* release memory                                                                              */
   for (i = 0; i < k; i++) {free(c1[i]);}
   free(c1);
   free(labels);
   free(cluster_score);
   free(cluster_id);
   free(tmp_center);

   printf("\tClustering done.\n");
   return 0;
}
