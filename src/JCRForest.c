#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include "utils.h"
#include <Rinternals.h>

double H_c(double *freqs, int *nclass){
  
  double score=0;
  
  for(int i = 0; i < *nclass; i++){
    if(freqs[i] > 0){
      score += freqs[i] * log(freqs[i]);
    }
  }
  return(score);
}

void find_best_split(double *x, double *yc, int *yf,int *nclass, int *mtry,int *nsample,int *nvar,int start, int end, int *ndind) {

  int var_ind[*nvar],j;
  
  for(int i=0; i<*nvar;i++) var_ind[i] = i;
  
  double best_score = 0.0;
  double curr_score = 0.0;
  double parent_score = 0.0;
  
  double pcx[*nclass];
  double pcxl[*nclass];
  double pcxr[*nclass];
  
  // adapt:: take only values from the given nodes
  for(int i=0; i < *nclass; i++) pcx[i] = 0;
  for(int i=start; i < end; i++) pcx[yf[ndind[i]]]++;
  for(int i=0; i < *nclass; i++) pcx[i]/= *nsample;
  
  int last = *nvar;
  double x_sort[*nsample];
  int x_sort_ind[*nsample];
  for(int i=0; i<*nsample; i++) x_sort_ind[i] = i;
  double yc_sorted[*nsample],yf_sorted[*nsample];
  
  for(int i=0; i < *mtry; i++){
    // select variable
    j = (int) floor(unif_rand() * (double) last);
    swap_integers(&var_ind[j],&var_ind[last-1]);
    last--;

    // copy the variables in temporary vectors
    for(int k=0; k< *nsample; k++) x_sort[k] = x[var_ind[last] * *nsample+k];  
    // sort x variables and obtain the ordering -- adapt x_sort_ind properly
    R_qsort_I(x_sort,x_sort_ind,1,*nsample);
    //printf("%d\n\n",var_ind[last]);
    // assign parent score
    
    // assign children distribution vectors
    memcpy(pcxl,pcx,*nclass);
    for(int k=0; k < *nclass; k++) pcxr[k] = 0;
    
    for(int k=start; k<end;k++){ // for each possible split
      //curr_score = H_c();
      if(curr_score-parent_score < best_score){
        best_score = curr_score - parent_score;
        
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