/*  
 *  Source code for Python module to be used in coexpression analysis 
 *
 *  Alan Gray, EPCC, The University of Edinburgh
 *  June 2013
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* forward declaration */
static void myRankData(double *data, double *sortedRanks, 
		       int *originalIndex, int size); 


/* main routine for calculating Pearson or Spearman coefficients */
void getPearson(double *dataAll, double *resultAll, int* idxa_all,
		int* idxb_all, int nPromPairs, int n){


  int ompNumThreads;
  double t1,t2,t3,t4;

  if(!getenv ("OMP_NUM_THREADS")) 
    {
      omp_set_num_threads(4);
      //printf("set the OMP_NUM_THREADS envronment variable to match the number of cores in your system \n");
    }

  
#pragma omp parallel
  {
    ompNumThreads=omp_get_num_threads();
  }
  
  //printf("Running getPearson function with %d OpenMP threads\n\n",ompNumThreads);

  // allocate memory
  double *adata= (double*) malloc(ompNumThreads*n*sizeof(double));
  double *bdata= (double*) malloc(ompNumThreads*n*sizeof(double)); 
  double *dworkspace= (double*) malloc(ompNumThreads*n*sizeof(double)); 
  int *iworkspace= (int*) malloc(ompNumThreads*n*sizeof(int)); 

  int x;

  #pragma omp parallel for 
  for(x=0; x<nPromPairs; x++){

    
    int threadId=omp_get_thread_num();
    
    int i;
        
    double xysum=0.;
    double xxsum=0.;
    double yysum=0.;    
    double adataMean=0.;
    double bdataMean=0.;
    double adatasum=0.;
    double bdatasum=0.;

    
    /* get promoter pair */
    int proma_idx=idxa_all[x];
    int promb_idx=idxb_all[x];

    //printf("\n");
   
    for (i=0; i<n; i++){
      adata[threadId*n+i]=dataAll[proma_idx*n+i];
      bdata[threadId*n+i]=dataAll[promb_idx*n+i];

      //   if (x==0){
      //printf("%f, ",adata[threadId*n+i]);
      //printf("%f, ",bdata[threadId*n+i]);
      //}

      adatasum+=adata[threadId*n+i];
      bdatasum+=bdata[threadId*n+i];    
    }
    //printf("\n");
    
    // if sum of data is zero, set result to NAN and skip this iteration 
    if (adatasum==0. || bdatasum==0.){
      resultAll[x]=NAN;
      continue;
    }
    
    // calculate means
    for (i=0; i<n; i++){
      adataMean+=adata[threadId*n+i];
      bdataMean+=bdata[threadId*n+i];
    }
    
    adataMean/=n;
    bdataMean/=n;
    
    
    // calculate pearson coefficient
    for (i=0; i<n; i++){
      xysum+=(adata[threadId*n+i]-adataMean)
	*(bdata[threadId*n+i]-bdataMean);
      xxsum+=pow(adata[threadId*n+i]-adataMean,2);
      yysum+=pow(bdata[threadId*n+i]-bdataMean,2);
    }
    
    // save result
    resultAll[x]=xysum/(sqrt(xxsum)*sqrt(yysum));
    
  }    
  
  
  // clean up 
  free(adata);
  free(bdata);
  free(dworkspace);
  free(iworkspace);

  //printf("%1.14e s full C\n", omp_get_wtime()-t1); t1=omp_get_wtime();
  
  return;
  

}