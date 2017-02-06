#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include "utils.h"

double gini_impurity(int nind, int nclass){
  
  double freq[nclass];
  for(int i=0; i<nclass; i++) freq[i] = 0;
  
  for(int i = 0; i < nind; i++){
    
  }
}

void find_best_split(double *x, double *yc, int *yf,int nclass, int mtry,int nsample,int nvar) {
  
  int var_ind[nvar],j;
  for(int i=0; i<nvar;i++) var_ind[i] = i;
  
  double best_score = 0.0;
  double curr_score = 0.0;
  double parent_score = 0.0;
  
  int last = nvar;
  double x_sort[nsample];
  int x_sort_ind[nsample];
  double yc_sorted[nsample],yf_sorted[nsample];
  
  for(int i=0; i < mtry; i++){
    // select variable
    j = (int) floor(unif_rand() * (last));
    swap_integers(&var_ind[j],&var_ind[last-1]);
    last--;
    for(int k=0; k<nsample; k++) x_sort[k] = x[var_ind[last] * nsample+k];  
    // sort x variables and obtain the ordering
    R_qsort_I(x_sort,x_sort_ind,1,nsample);
    
    // assign parent score
    
    for(int k=0; k<(nsample-1);k++){ // for each possible split
    
      if(curr_score-parent_score < best_score){
        best_score = curr_score - parent_score;
        
      }  
    }
    // the following is not completely necessary, one can also take these values out on the go
    // [however, that works ]
    /*for(int k=0; k<nsample; k++){
      yc_sorted[k] = yc[x_sort_ind[k]];
      yf_sorted[k] = yf[x_sort_ind[k]];
    }*/
  }
  
  
  
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int nclass, int curr_tree, int *ntree, int *nrnodes,int *ldaughter, int *rdaughter,
                    int mtry,int nsample,int nvar) {
  
  
  for(int i=0; i < *nrnodes; i++){
    find_best_split(x,yc,yf,mtry,nsample,nvar);
  }
  // create vector containing which node is whose left and whose right daughter
  ldaughter[curr_tree] = curr_tree;
  rdaughter[0] = 11;
  
}


void build_jcr_forest(double *x, double* yc, int* yf, int *nclass, int *nsample , int *nvar, int *mtry, int *ntree, 
                      int *nrnodes, int *ldaughter, int *rdaughter, int *node_status, int *node_var,
                      int *node_xvar) {
  

  GetRNGstate();
  //printf("value is %i",*ntree);
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
    //printf("the random number %d \n",ind);
    // actual tree building
    build_jcr_tree(x_bag,yc_bag,yf_bag,*nclass,i,ntree,nrnodes,ldaughter+idx,rdaughter+idx,*mtry,*nsample,*nvar);
    
    // tree prediction
  }
  
  PutRNGstate();  
  
}