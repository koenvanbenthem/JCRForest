// code

void forestpred(double *x, double *yc, int *yf, int *nclass, int *nsample, int *ldaughter, int *rdaughter, int *yf_pred, int *ntree, int *nrnodes, 
                int *node_var, double *node_xvar){
    

    for(int i=0; i < *nsample; i++){
      
      int yf_p[*nclass];
      for(int k=0; k < *nclass; k++){ yf_p[k]=0; }
      
      for(int j=0; j < *ntree; j++){
        
        int found = 0;
        int curr_ind = 0;
        
        while(found == 0){
          
          if(x[i + *nsample * node_var[curr_ind + *nrnodes * j]] > node_xvar[curr_ind + *nrnodes * j]){
            curr_ind = ldaughter[curr_ind + j * *nrnodes];
            //printf("left%d\t\t",curr_ind);
          }else{
            curr_ind = rdaughter[curr_ind + j * *nrnodes];
            //printf("right%d\t\t",curr_ind);
          }
          //found=1;
          printf("%d\t%d\t%d\t%d\n",curr_ind,j,*nrnodes,(curr_ind + j * *nrnodes));
          if(ldaughter[curr_ind + j * *nrnodes] == 0){
            found = 1;
            yf[*nclass * i + yf_pred[curr_ind + j * *nrnodes]]++;
            printf("in the loop\n");
          }
          
        }
      }
    }
}