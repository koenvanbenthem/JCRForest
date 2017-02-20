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

void find_best_split(double *x, double *yc, int *yf,int *nclass, int *mtry,int *nsample,int *nvar,
                     int start, int end, int *ndind, int *best_var, double *best_split, int *best_k, int *yf_predr, int *yf_predl, 
                     double *yc_mu_predr, double *yc_mu_predl, double *yc_sd_predr, double *yc_sd_predl, FILE *fp) {
  
  // vector for variable indices - variables of choice are drawn from this vector
  int var_ind[*nvar];
  for(int i=0; i<*nvar;i++) var_ind[i] = i;
  int last = *nvar;
  *best_var = -1;
  *best_split = -1.0;
  *best_k = -1;
  
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

  parent_score = H_c(pcx,nclass,end-start);
  
  for(int i=0; i < *mtry; i++){

    // select variable
    int j = (int) floor(unif_rand() * (double) last);
    swap_integers(&var_ind[j],&var_ind[last-1]);
    last--;

    // copy the variables in temporary vectors
    for(int k=start; k< end; k++) x_sort[k] = x[var_ind[last] * *nsample+ndind[k]];  
    for(int k=start; k<end; k++) x_sort_ind[k] = ndind[k];
    
    // sort x variables and obtain the ordering -- adapt x_sort_ind properly
    fprintf(fp,"[%d,%d]\t%d \t %f %f ;;; trying var %d\n",start+1,end,*nsample,x_sort[start],x_sort[end-1],var_ind[last]);
    fflush(fp);
    R_qsort_I(x_sort,x_sort_ind,start+1,end);

    // assign parent score
    
    // assign children distribution vectors (left child gets all, right gets none)
    memcpy(pcxl,pcx,*nclass * sizeof(int));
    
    //printf("%d, %d\n",pcxl[0],pcxl[1]);
    for(int k=0; k < *nclass; k++) pcxr[k] = 0;

    //printf("%f\n",parent_score);
    //curr_score = parent_score; // initially one of the children is the same as the parent, the other is empty
    
    int Nl = end-start;
    int Nr = 0;
    for(int k=start; k<(end-1); k++){ // for each possible split
      
      pcxl[yf[x_sort_ind[k]]-1]--;
      pcxr[yf[x_sort_ind[k]]-1]++;
      Nr++;
      Nl--;
      curr_score = ((double) Nr/((double) (end-start)))*H_c(pcxr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_c(pcxl,nclass,Nl);
      double nvv = curr_score;
      /*printf("(%d, %d)\t",pcxl[0],pcxl[1]);
      printf("(%d, %d)\t",pcxr[0],pcxr[1]);
      printf("%f \n",curr_score);*/
      //printf("pcxl1 = %d, pcxl2 = %d, pcxr1 = %d, pcxr2 = %d, Nr = %d, Nl = %d, k = %d, var = %d, the newest score is %f\n",pcxl[0],pcxl[1],pcxr[0],pcxr[1],Nr,Nl,k,var_ind[last],parent_score - nvv);
      if(parent_score-curr_score > best_score){
        best_score = parent_score-curr_score;
        *best_var = var_ind[last];
        *best_k = k;
        //printf("Bla");
        
        *yc_mu_predr = 9.5;
        *yc_mu_predl = 1.8;
        
        *yc_sd_predr = 1.25;
        *yc_sd_predl = 1.12; 
        
        *yf_predl = which_max(pcxl,*nclass);
        *yf_predr = which_max(pcxr,*nclass);
        
      }  
    }

  }
  
  if(*best_k != -1){
    // store variable
    for(int k=start; k < end; k++){
      x_sort[k] = x[*best_var * *nsample+ndind[k]];  
      //fprintf(fp,"var %d, k %d index %d\n",*best_var,k,ndind[k]);
    }
    // storing the new order of ndind
    fprintf(fp,"in if [%d,%d]\t%d\t%f %f (bestk) %d\n",start+1,end,*nsample,x_sort[start],x_sort[end-1],*best_k);
    fflush(fp);
    R_qsort_I(x_sort,ndind,start+1,end);
    /*if(*best_k >= end){
      printf("Quoi");
    }*/
    *best_split = 0.5*x_sort[*best_k]+0.5*x_sort[*best_k+1];
  }
  //printf("[%d,%d] Best variable is... %d with a score of %f (ps=%f)\n",start,end,best_var,best_score,parent_score);
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int *nclass, int curr_tree, int *ntree, int *nrnodes, int *minsize, int *ldaughter, int *rdaughter,
                    int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_var, double *node_xvar, int *mtry,int *nsample,int *nvar, FILE *fp) {
  
  int ndstart[*nrnodes];
  int ndend[*nrnodes];
  int ndind[*nsample];
  int last_node = 0;
  ndstart[0] = 0;
  ndend[0] = *nsample;
  for(int i=0; i<*nsample; i++) ndind[i] = i;
  
  
  int best_var, best_k,yf_predr,yf_predl;
  double best_split,yc_mu_predr,yc_mu_predl,yc_sd_predr,yc_sd_predl;
  
  for(int i=0; i < *nrnodes; i++){//*nrnodes; i++){
    
    if(last_node > *nrnodes-3 || i > last_node) break;
    
    if(ndend[i] - ndstart[i] <= (*minsize - 1)){
      //printf("[%d,%d] nothing left to do \n",ndstart[i],ndend[i]);
      continue;
    }
    yf_predr = -2;
    yf_predl = -2;
    //int i=0; // temporary for testing purposes
    find_best_split(x,yc,yf,nclass,mtry,nsample,nvar,ndstart[i],ndend[i],ndind, &best_var, &best_split, &best_k, &yf_predr, &yf_predl,&yc_mu_predr,&yc_mu_predl,&yc_sd_predr,&yc_sd_predl,fp);
    
    node_var[i] = best_var;
    node_xvar[i] = best_split;
    // saving these best entries
    if(best_split != -1 & best_var != -1.0){
      //printf("[%d,%d]\t Setting node_var[%d] to %d \t\t with a value of %f\t\t best k is %d\n",ndstart[i],ndend[i],i,best_var,best_split,best_k);
      //printf("[%d,%d] --> [%d,%d] u [%d,%d]",)
      ndstart[last_node+1] = ndstart[i];
      ndend[last_node+1] = best_k+1;
      ndstart[last_node+2] = best_k+1;
      ndend[last_node+2] = ndend[i]; 
      rdaughter[i] = last_node + 1;
      ldaughter[i] = last_node + 2;
      
      yf_pred[last_node+2] = yf_predl;
      yf_pred[last_node+1] = yf_predr;
      
      yc_mu_pred[last_node+2] = yc_mu_predl;
      yc_mu_pred[last_node+1] = yc_mu_predr;
      
      yc_sd_pred[last_node+2] = yc_sd_predl;
      yc_sd_pred[last_node+1] = yc_sd_predr;
      
      last_node = last_node + 2;
    }else{
      //printf("gnagna\n");
    }
    
  }
  

}

void build_jcr_forest(double *x, double* yc, int* yf, int *nclass, int *nsample , int *nvar, int *mtry, int *ntree, 
                      int *nrnodes, int *minsize, int *ldaughter, int *rdaughter, int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_status, 
                      int *node_var, double *node_xvar, double *dum_vect, int *dum_long, int *dum_ind) {
  

  GetRNGstate();
  
  FILE *fp;
  fp = fopen("log.txt","w");
  //fprintf(fp,"blabla %d",7);

  double ran_num;
  
  for(int i=0; i < *nrnodes; i++){
    for(int j=0; j < *ntree; j++){
      node_var[j * *nrnodes + i] = -1;
      yf_pred[j * *nrnodes + i] = -1;
    }
  }
  double x_bag[*nsample * *nvar];
  double yc_bag[*nsample];
  int yf_bag[*nsample];
  int ind;
  // Building the trees
  for(int i = 0; i < *ntree; i++){
    printf("A NEW TREE IS BORN\n-----------------------------------\n");
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
    build_jcr_tree(x_bag,yc_bag,yf_bag,nclass,i,ntree,nrnodes,minsize,ldaughter+idx,rdaughter+idx,yf_pred+idx,yc_mu_pred+idx,yc_sd_pred+idx,node_var+idx,node_xvar+idx,mtry,nsample,nvar,fp);
    
    // tree prediction
  }

  fclose(fp);
  PutRNGstate(); 
  
}