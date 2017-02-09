#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include "utils.h"
#include <Rinternals.h>

double H_c(int *counts, int *nclass, int N){
  
  double score=0;
  
  for(int i = 0; i < *nclass; i++){
    if(counts[i] > 0){
      score -= (double) counts[i]/ (double) N * log((double) counts[i]/ (double) N);
    }
  }
  return(score);
}

void find_best_split(double *x, double *yc, int *yf,int *nclass, int *mtry,int *nsample,int *nvar,int start, int end, int *ndind) {

  // vector for variable indices - variables of choice are drawn from this vector
  int var_ind[*nvar];
  for(int i=0; i<*nvar;i++) var_ind[i] = i;
  int last = *nvar;
  
  // scoring variables
  double best_score = 0.0;
  double curr_score = 0.0;
  double parent_score = 0.0;
  
  // keeping track of the frequencies of different classes
  // parent, left child, right child
  int pcx[*nclass];
  int pcxl[*nclass];
  int pcxr[*nclass];
  
   // setting initial values (i.e. parent frequencies)
  for(int i=0; i < *nclass; i++) pcx[i] = 0;
  for(int i=start; i < end; i++) pcx[yf[ndind[i]]-1]++;
  //for(int i=0; i < *nclass; i++) pcx[i]/= (end-start);

  // vectors for storing the explanatory variable (x_choice) as well as the associated indices
  double x_sort[*nsample];
  int x_sort_ind[*nsample];
  for(int i=start; i<end; i++) x_sort_ind[i] = ndind[i];
  
  double yc_sorted[*nsample];
  int yf_sorted[*nsample];

  for(int i=0; i < *mtry; i++){

    // select variable
    int j = (int) floor(unif_rand() * (double) last);
    swap_integers(&var_ind[j],&var_ind[last-1]);
    last--;

    // copy the variables in temporary vectors
    for(int k=start; k< end; k++) x_sort[k] = x[var_ind[last] * *nsample+k];  
    // sort x variables and obtain the ordering -- adapt x_sort_ind properly
    R_qsort_I(x_sort,x_sort_ind,start+1,end);

    // assign parent score
    
    // assign children distribution vectors (left child gets all, right gets none)
    memcpy(pcxl,pcx,*nclass * sizeof(double));
    //printf("%d, %d\n",pcxl[0],pcxl[1]);
    for(int k=0; k < *nclass; k++) pcxr[k] = 0;
    parent_score = H_c(pcx,nclass,end-start);
    printf("%f\n",parent_score);
    //curr_score = parent_score; // initially one of the children is the same as the parent, the other is empty
    
    int Nl = end-start;
    int Nr = 0;
    for(int k=start; k<(end-1);k++){ // for each possible split
      
      pcxl[yf[x_sort_ind[k]]-1]--;
      pcxr[yf[x_sort_ind[k]]-1]++;
      Nr++;
      Nl--;
      curr_score = ((double) Nr/((double) (end-start)))*H_c(pcxr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_c(pcxl,nclass,Nl);
      double nvv = curr_score;
      printf("pcxl1 = %d, pcxl2 = %d, pcxr1 = %d, pcxr2 = %d, Nr = %d, Nl = %d, k = %d, var = %d, the newest score is %f\n",pcxl[0],pcxl[1],pcxr[0],pcxr[1],Nr,Nl,k,var_ind[last],nvv);
      if(parent_score-curr_score < best_score){
        best_score = parent_score-curr_score;
        
      }  
    }

  }
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int *nclass, int curr_tree, int *ntree, int *nrnodes,int *ldaughter, int *rdaughter,
                    int *mtry,int *nsample,int *nvar) {
  
  int ndstart[*nrnodes];
  int ndend[*nrnodes];
  int ndind[*nsample];
  ndstart[0] = 0;
  ndend[0] = *nsample;
  for(int i=0; i<*nsample; i++) ndind[i] = i;
  
  //for(int i=0; i < *nrnodes; i++){
    int i=0; // temporary for testing purposes
    find_best_split(x,yc,yf,nclass,mtry,nsample,nvar,ndstart[i],ndend[i],ndind);
  //}
  
  // create vector containing which node is whose left and whose right daughter
  ldaughter[curr_tree] = curr_tree;
  rdaughter[0] = 11;
  
}

void build_jcr_forest(double *x, double* yc, int* yf, int *nclass, int *nsample , int *nvar, int *mtry, int *ntree, 
                      int *nrnodes, int *ldaughter, int *rdaughter, int *node_status, int *node_var,
                      int *node_xvar, double *dum_vect, int *dum_long, int *dum_ind) {
  

  GetRNGstate();
  
  double ran_num;
  
  double x_bag[*nsample * *nvar];
  double yc_bag[*nsample];
  int yf_bag[*nsample];
  int ind;
  // Building the trees
  for(int i = 0; i < *ntree; i++){
    
    int idx = i * *nrnodes; // coordinates for the specific tree 
    
    // data selection - only with replacement
    for(int j=0; j < *nsample; j++){
      ran_num = unif_rand();
      ind = (int) floor(ran_num * *nsample);
      yf_bag[j] = yf[ind];
      yc_bag[j] = yc[ind];
      
      for(int k=0; k < *nvar; k++){
        x_bag[k * *nsample + j] = x[k * *nsample + ind];
      }
    }
    
    // actual tree building
    build_jcr_tree(x_bag,yc_bag,yf_bag,nclass,i,ntree,nrnodes,ldaughter+idx,rdaughter+idx,mtry,nsample,nvar);
    
    // tree prediction
  }

  R_qsort_I(dum_vect,dum_ind,5,10);

  PutRNGstate(); 
  
}