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

double H_rc(int *counts,double *sd, int *nclass, int N){
  double score = 0;
  
  for(int i=0; i < *nclass; i++){
    if(counts[i] > 0 & sd[i] > 0){
      score += ((double) counts[i]) * 0.5*log(2.0*PI*exp(1)*sd[i])/((double) N); 
    }
  }
  return score;
}

void robust_mean(int *N, double kappa, int nclass, double *mean_child, double *mean_parent, double *store){
  for(int i=0; i<nclass; i++){
    store[i] = ((double) N[i] / (kappa + (double) N[i])) * mean_child[i] + (kappa/(kappa + (double) N[i])) * mean_parent[i];
  }
}

void robust_sd(int *N, double kappa, double nu, double n, int nclass, double *sd_child, double *sd_parent,double *store, double *mean_child, double *mean_parent){
  for(int i=0; i<nclass; i++){
    double Z = nu + n - 1 + N[i];
    store[i] = ((double) N[i] / Z) * sd_child[i] + (nu+n-1)/Z * sd_parent[i] + (kappa * (double) N[i]) / (Z * (kappa + N[i])) * pow(mean_parent[i] - mean_child[i],2.0);
  }
}

void find_best_split(double *x, double *yc, int *yf,int *nclass, int *mtry,int *nsample,int *nvar, int *minsize, double H0_c, double H0_rc,
                     int start, int end, int *ndind, int *best_var, double *best_split, int *best_k, int *yf_predr, int *yf_predl, 
                     double *yc_mu_predr, double *yc_mu_predl, double *yc_sd_predr, double *yc_sd_predl, double *kappa, double *nu,
                     double *y_mu_c,double *y_sd_c,int pid, double *best_meanl, double *best_meanr, double *best_sdl, double *best_sdr) {
  

  // vector for variable indices - variables of choice are drawn from this vector
  int var_ind[*nvar];
  for(int i=0; i<*nvar;i++) var_ind[i] = i;
  int last = *nvar;
  *best_var = -1;
  *best_split = -1.0;
  *best_k = -1;
  
  // scoring variables
  double best_score = 0.0;
  double curr_score_disc = 0.0;
  double curr_score_cont = 0.0;
  double parent_score = 0.0;
  
  // keeping track of the frequencies of different classes
  // parent, left child, right child
  int pcx[*nclass];
  int pcxl[*nclass];
  int pcxr[*nclass];
  
  // Continuous variables
  double meanl[*nclass];
  double meanr[*nclass];
  double sdl[*nclass];
  double sdr[*nclass];
  
  double real_meanl[*nclass];
  double real_meanr[*nclass];
  double real_sdl[*nclass];
  double real_sdr[*nclass];
  
   // setting initial values (i.e. parent frequencies)
  for(int i=0; i < *nclass; i++) pcx[i] = 0;
  for(int i=start; i < end; i++) pcx[yf[ndind[i]]-1]++;
  //for(int i=0; i < *nclass; i++) pcx[i]/= (end-start);

  // vectors for storing the explanatory variable (x_choice) as well as the associated indices
  double x_sort[*nsample];
  int x_sort_ind[*nsample];
  for(int i=start; i<end; i++) x_sort_ind[i] = ndind[i];
  
  //double yc_sorted[*nsample];
  //int yf_sorted[*nsample];
  if(pid==-1){
    calc_mean_class_I(yc,yf,nclass,start,end,x_sort_ind,meanr);
    calc_sd_class_I(yc,yf,nclass,start,end,x_sort_ind,sdr,meanr);
    parent_score = 0.5 * (H_c(pcx,nclass,end-start)/H0_c + H_rc(pcx,sdr,nclass,end-start)/H0_rc);
  } else {
    parent_score = 0.5 * (H_c(pcx,nclass,end-start)/H0_c + H_rc(pcx,y_sd_c+pid*(*nclass),nclass,end-start)/H0_rc);
  }
  pid = -1;
  for(int i=0; i < *mtry; i++){

    // select variable
    int j = (int) floor(unif_rand() * (double) last);
    swap_integers(&var_ind[j],&var_ind[last-1]);
    last--;

    // copy the variables in temporary vectors
    for(int k=start; k< end; k++) x_sort[k] = x[var_ind[last] * *nsample+ndind[k]];  
    for(int k=start; k<end; k++) x_sort_ind[k] = ndind[k];
    
    // sort x variables and obtain the ordering -- adapt x_sort_ind properly
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
    
    for(int k=start; k < start+*minsize; k++){
      pcxl[yf[x_sort_ind[k]]-1]--;
      pcxr[yf[x_sort_ind[k]]-1]++;
      Nr++;
      Nl--;  
    }
    
    for(int k=start+*minsize; k<(end-1-*minsize); k++){ // for each possible split
      
      pcxl[yf[x_sort_ind[k]]-1]--;
      pcxr[yf[x_sort_ind[k]]-1]++;
      Nr++;
      Nl--;
      
      // make robust
      calc_mean_class_I(yc,yf,nclass,start,k+1,x_sort_ind,meanr);
      calc_sd_class_I(yc,yf,nclass,start,k+1,x_sort_ind,sdr,meanr);
      
      calc_mean_class_I(yc,yf,nclass,k+1,end,x_sort_ind,meanl);
      calc_sd_class_I(yc,yf,nclass,k+1,end,x_sort_ind,sdl,meanl);
      double n = 1;
      if(pid!=-1){
        
        robust_mean(pcxl,*kappa,*nclass,meanl,y_mu_c+(pid* *nclass),real_meanl);
        robust_mean(pcxr,*kappa,*nclass,meanr,y_mu_c+(pid* *nclass),real_meanr);
        
        robust_sd(pcxl,*kappa,*nu,n,*nclass,sdl,y_sd_c+(pid * *nclass),real_sdl,meanl,y_mu_c+(pid* *nclass));
        robust_sd(pcxr,*kappa,*nu,n,*nclass,sdr,y_sd_c+(pid * *nclass),real_sdr,meanr,y_mu_c+(pid* *nclass));
      } else {
        memcpy(real_sdr,sdr,*nclass * sizeof(double));
        memcpy(real_sdl,sdl,*nclass * sizeof(double));
        memcpy(real_meanl,meanl,*nclass * sizeof(double));
        memcpy(real_meanr,meanr,*nclass * sizeof(double));
      }
      curr_score_cont =  ((double) Nr/((double) (end-start)))*H_rc(pcxr,real_sdr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_rc(pcxl,real_sdl,nclass,Nl);
      curr_score_disc = ((double) Nr/((double) (end-start)))*H_c(pcxr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_c(pcxl,nclass,Nl);
      double nvv = 0.5 * (curr_score_disc/H0_c + curr_score_cont/H0_rc);

      if(parent_score-nvv > best_score){
        best_score = parent_score-nvv;
        *best_var = var_ind[last];
        *best_k = k;

        // Assign parent robust mean and sd here
        memcpy(best_meanl,real_meanl,*nclass * sizeof(double));
        memcpy(best_meanr,real_meanr,*nclass * sizeof(double));
        memcpy(best_sdl,real_sdl,*nclass * sizeof(double));
        memcpy(best_sdr,real_sdr,*nclass * sizeof(double));
        
        
        *yc_mu_predr = weighted_average(real_meanr,pcxr,*nclass);
        *yc_mu_predl = weighted_average(real_meanl,pcxl,*nclass);
        
        *yc_sd_predr = 1.25;
        *yc_sd_predl = 1.12; 
        
        *yf_predl = which_max(pcxl,*nclass)+1;
        *yf_predr = which_max(pcxr,*nclass)+1;
      }  
    }

  }
  
  if(*best_k != -1){
    // store variable
    for(int k=start; k < end; k++){
      x_sort[k] = x[*best_var * *nsample+ndind[k]];  
    }
    // storing the new order of ndind
    R_qsort_I(x_sort,ndind,start+1,end);

    *best_split = 0.5*x_sort[*best_k]+0.5*x_sort[*best_k+1];
  }
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int *nclass, int curr_tree, int *ntree, int *nrnodes, int *minsize, int *ldaughter, int *rdaughter,
                    int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_var, double *node_xvar, int *mtry,int *nsample,int *nvar, double *kappa, double *nu) {
  
  int ndstart[*nrnodes];
  int ndend[*nrnodes];
  int ndind[*nsample];
  int last_node = 0;
  int parent_id[*nrnodes];
  for(int i=0; i<*nrnodes;i++) parent_id[i] = -1;
  ndstart[0] = 0;
  ndend[0] = *nsample;
  for(int i=0; i<*nsample; i++) ndind[i] = i;
  
  int best_var, best_k,yf_predr,yf_predl;
  double best_split,yc_mu_predr,yc_mu_predl,yc_sd_predr,yc_sd_predl;
  
  double y_mu_c[*nclass * *nrnodes];
  double y_sd_c[*nclass * *nrnodes];
  // void calc_sd_class(double *vector, int *class_vector, int *nclass, int nsample, double *store, double *mean){
  //void calc_mean_class(double *vector, int *class_vector, int *nclass, int nsample,double *store)
  
  calc_mean_class(yc,yf,nclass,*nsample,y_mu_c);
  calc_sd_class(yc,yf,nclass,*nsample,y_sd_c,y_mu_c);

  int pcx[*nclass];
  for(int i=0; i < *nclass; i++) pcx[i] = 0;
  for(int i=0; i < *nsample; i++) pcx[yf[ndind[i]]-1]++;
  
  yf_pred[0] = which_max(pcx,*nclass)+1;
  yc_mu_pred[0] = weighted_average(y_mu_c,pcx,*nclass);
  
  double H0_c = H_c(pcx,nclass,*nsample);
  double H0_rc = H_rc(pcx,y_sd_c,nclass,*nsample);
  
  double best_meanl[*nclass];
  double best_meanr[*nclass];
  double best_sdl[*nclass];
  double best_sdr[*nclass];
  
  for(int i=0; i < *nrnodes; i++){//*nrnodes; i++){
    
    if(last_node > *nrnodes-3 || i > last_node) break;
    
    if(ndend[i] - ndstart[i] <= (*minsize * 2 - 1)){
      continue;
    }
    yf_predr = -2;
    yf_predl = -2;
    //int i=0; // temporary for testing purposes
    find_best_split(x,yc,yf,nclass,mtry,nsample,nvar,minsize,H0_c,H0_rc,ndstart[i],ndend[i],ndind, &best_var, 
                    &best_split, &best_k, &yf_predr, &yf_predl,&yc_mu_predr,&yc_mu_predl,
                    &yc_sd_predr,&yc_sd_predl,kappa,nu,y_mu_c,y_sd_c,parent_id[i],best_meanl,best_meanr,best_sdl,best_sdr);
    
    node_var[i] = best_var;
    node_xvar[i] = best_split;


    
    // saving these best entries
    if(best_split != -1 & best_var != -1.0){
      
      ndstart[last_node+1] = ndstart[i];
      ndend[last_node+1] = best_k+1;
      ndstart[last_node+2] = best_k+1;
      ndend[last_node+2] = ndend[i]; 
      rdaughter[i] = last_node + 1;
      ldaughter[i] = last_node + 2;
      
      memcpy(y_mu_c+((last_node+1) * *nclass),best_meanr,*nclass * sizeof(double));
      memcpy(y_sd_c+((last_node+1) * *nclass),best_sdr,*nclass * sizeof(double));
      memcpy(y_mu_c+((last_node+2) * *nclass),best_meanl,*nclass * sizeof(double));
      memcpy(y_sd_c+((last_node+2) * *nclass),best_sdl,*nclass * sizeof(double));
      
      parent_id[last_node+1] = i;
      parent_id[last_node+2] = i;
      
      yf_pred[last_node+2] = yf_predl;
      yf_pred[last_node+1] = yf_predr;
      
      yc_mu_pred[last_node+2] = yc_mu_predl;
      yc_mu_pred[last_node+1] = yc_mu_predr;
      
      yc_sd_pred[last_node+2] = yc_sd_predl;
      yc_sd_pred[last_node+1] = yc_sd_predr;
      
      last_node = last_node + 2;
    }
  }
}

void build_jcr_forest(double *x, double* yc, int* yf, int *nclass, int *nsample , int *nvar, int *mtry, int *ntree, 
                      int *nrnodes, int *minsize, int *ldaughter, int *rdaughter, int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_status, 
                      int *node_var, double *node_xvar, double *dum_vect, int *dum_long, int *dum_ind, double *kappa,double *nu) {
  

  GetRNGstate();
  
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
    build_jcr_tree(x_bag,yc_bag,yf_bag,nclass,i,ntree,nrnodes,minsize,ldaughter+idx,rdaughter+idx,yf_pred+idx,yc_mu_pred+idx,yc_sd_pred+idx,node_var+idx,node_xvar+idx,mtry,nsample,nvar,kappa,nu);
    
  }
  
  PutRNGstate(); 
  
}