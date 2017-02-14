// code

void forestpred(double *x, double *yc, int *yf, int *nsample, int *ldaughter, int *rdaughter, int *ntree, int *nrnodes, 
                int *node_var, double *node_xvar){
    
    
    for(int i=0; i < *nsample; i++){
      for(int j=0; j < *ntree; j++){
        
        int found = 0;
        int curr_ind = 0;
        
        while(found == 0){
          
          if(x[i + *nsample * node_var[curr_ind + *ntree * j]] < node_xvar[curr_ind + *ntree * j]){
            curr_ind = ldaughter[curr_ind + j * *nrnodes];
          }else{
            curr_ind = rdaughter[curr_ind + j * *nrnodes];
          }
          //found=1;
          printf("%d\t%d\t%d\n",curr_ind,j,*nrnodes);
          if(ldaughter[curr_ind + j * *nrnodes] == 0){
            found = 1;
            printf("in the loop\n");
          }
          
        }
      }
    }
}