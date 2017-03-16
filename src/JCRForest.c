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

void robust_mean(int N, double kappa, int nclass, double *mean_child, double *mean_parent, double *store){
  for(int i=0; i<nclass; i++){
    store[i] = ((double) N / (kappa + (double) N)) * mean_child[i] + (kappa/(kappa + (double) N)) * mean_parent[i];
  }
}

void robust_sd(int N, double kappa, double nu, double n, int nclass, double *sd_child, double *sd_parent,double *store, double *mean_child, double *mean_parent){
  double Z = nu + n - 1 + N;

  for(int i=0; i<nclass; i++){
    store[i] = ((double) N / Z) * sd_child[i] + (nu+n-1)/Z * sd_parent[i] + (kappa * (double) N) / (Z * (kappa + N)) * pow(mean_parent - mean_child,2.0);
  }
}

void find_best_split(double *x, double *yc, int *yf,int *nclass, int *mtry,int *nsample,int *nvar, int *minsize, double H0_c, double H0_rc,
                     int start, int end, int *ndind, int *best_var, double *best_split, int *best_k, int *yf_predr, int *yf_predl, 
                     double *yc_mu_predr, double *yc_mu_predl, double *yc_sd_predr, double *yc_sd_predl, FILE *fp, double *kappa, double *nu,
                     double *y_mu_c,double *y_sd_c,int pid) {
  

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
      
      if(pid!=-1){
        //void robust_mean(int N, double kappa, int nclass, double *mean_child, double *mean_parent, double *store){
        robust_mean(Nl,*kappa,*nclass,meanl,y_mu_c+pid* *nclass,real_meanl);
        robust_mean(Nl,*kappa,*nclass,meanr,y_mu_c+pid* *nclass,real_meanr);
        robust_sd();
        robust_sd();
        //real_sdl = ;
        //real_sdr = ;
      }
      curr_score_cont =  ((double) Nr/((double) (end-start)))*H_rc(pcxr,sdr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_rc(pcxl,sdl,nclass,Nl);
      curr_score_disc = ((double) Nr/((double) (end-start)))*H_c(pcxr,nclass,Nr) + ((double) Nl/((double) (end-start)))*H_c(pcxl,nclass,Nl);
      //printf("%f, %f \n\n",curr_score_cont,curr_score_disc);
      double nvv = 0.5 * (curr_score_disc/H0_c + curr_score_cont/H0_rc);
      /*printf("(%d, %d)\t",pcxl[0],pcxl[1]);
      printf("(%d, %d)\t",pcxr[0],pcxr[1]);
      printf("%f \n",curr_score);*/
      //printf("pcxl1 = %d, pcxl2 = %d, pcxr1 = %d, pcxr2 = %d, Nr = %d, Nl = %d, k = %d, var = %d, the newest score is %f\n",pcxl[0],pcxl[1],pcxr[0],pcxr[1],Nr,Nl,k,var_ind[last],parent_score - nvv);
      if(parent_score-nvv > best_score){
        best_score = parent_score-nvv;
        *best_var = var_ind[last];
        *best_k = k;
        //printf("Bla");
        //print_array_int(pcxr,*nclass);
        //print_array_double(meanr,*nclass);
        *yc_mu_predr = weighted_average(meanr,pcxr,*nclass);
        //printf("Leads to prediction %f\n",*yc_mu_predr);
        *yc_mu_predl = weighted_average(meanl,pcxl,*nclass);
        
        *yc_sd_predr = 1.25;
        *yc_sd_predl = 1.12; 
        
        *yf_predl = which_max(pcxl,*nclass)+1;
        *yf_predr = which_max(pcxr,*nclass)+1;
        if(*yf_predl < 0 || *yf_predr < 0){
          printf("o,o");
        }
        
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

    *best_split = 0.5*x_sort[*best_k]+0.5*x_sort[*best_k+1];
  }
  //printf("[%d,%d] Best variable is... %d with a score of %f (ps=%f)\n",start,end,best_var,best_score,parent_score);
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int *nclass, int curr_tree, int *ntree, int *nrnodes, int *minsize, int *ldaughter, int *rdaughter,
                    int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_var, double *node_xvar, int *mtry,int *nsample,int *nvar, FILE *fp,double *kappa, double *nu) {
  
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
  
  double H0_c = H_c(pcx,nclass,*nsample);
  double H0_rc = H_rc(pcx,y_sd_c,nclass,*nsample);
  
  for(int i=0; i < *nrnodes; i++){//*nrnodes; i++){
    
    if(last_node > *nrnodes-3 || i > last_node) break;
    
    if(ndend[i] - ndstart[i] <= (*minsize - 1)/2){
      //printf("[%d,%d] nothing left to do \n",ndstart[i],ndend[i]);
      continue;
    }
    yf_predr = -2;
    yf_predl = -2;
    //int i=0; // temporary for testing purposes
    find_best_split(x,yc,yf,nclass,mtry,nsample,nvar,minsize,H0_c,H0_rc,ndstart[i],ndend[i],ndind, &best_var, 
                    &best_split, &best_k, &yf_predr, &yf_predl,&yc_mu_predr,&yc_mu_predl,
                    &yc_sd_predr,&yc_sd_predl,fp,kappa,nu,y_mu_c,y_sd_c,parent_id[i]);
    
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
      
      parent_id[last_node+1] = i;
      parent_id[last_node+2] = i;
      
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
  print_array_int(parent_id,*nrnodes);
  

}

void build_jcr_forest(double *x, double* yc, int* yf, int *nclass, int *nsample , int *nvar, int *mtry, int *ntree, 
                      int *nrnodes, int *minsize, int *ldaughter, int *rdaughter, int *yf_pred, double *yc_mu_pred, double *yc_sd_pred, int *node_status, 
                      int *node_var, double *node_xvar, double *dum_vect, int *dum_long, int *dum_ind, double *kappa,double *nu) {
  

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
      ind = j;//(int) floor(ran_num * *nsample);
      yf_bag[j] = yf[ind];
      yc_bag[j] = yc[ind];
      
      for(int k=0; k < *nvar; k++){
        x_bag[k * *nsample + j] = x[k * *nsample + ind];
      }
    }
    
    // actual tree building
    build_jcr_tree(x_bag,yc_bag,yf_bag,nclass,i,ntree,nrnodes,minsize,ldaughter+idx,rdaughter+idx,yf_pred+idx,yc_mu_pred+idx,yc_sd_pred+idx,node_var+idx,node_xvar+idx,mtry,nsample,nvar,fp,kappa,nu);
    
    // tree prediction
  }
  
  fclose(fp);
  PutRNGstate(); 
  
}