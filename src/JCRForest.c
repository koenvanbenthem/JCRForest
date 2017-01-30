#include <R.h>
#include <Rmath.h>
#include <stdio.h>

void calc_split() {
  
}

void build_jcr_tree(double *x, double *yc, int *yf, int curr_tree, int *ntree, int *nrnodes,int *ldaughter, int *rdaughter) {
  
  // create vector containing which node is whose left and whose right daughter
  ldaughter[curr_tree] = curr_tree;
  rdaughter[0] = 11;
  
}


void build_jcr_forest(double *x, double* yc, int* yf, int *nsample , int *nvar, int *ntree, int *nrnodes, int *ldaughter, 
                      int *rdaughter, int *node_status, int *node_var,int *node_xvar) {
  
  int xdim = *nsample;

  //printf("value is %i",*ntree);
  
  // Building the trees
  for(int i = 0; i < *ntree; i++){
    
    int idx = i * *nrnodes;
    
    // data selection
    
    // actual building
    build_jcr_tree(x,yc,yf,i,ntree,nrnodes,ldaughter+idx,rdaughter+idx);
  }
}
