#include <R.h>
#include <Rmath.h>

void build_jcr_tree(int* arr,int maxnodes) {
  
  // create vector containing which node is whose left and whose right daughter
  for(int i=0; i<maxnodes; i++){
    *(arr+i) = 2;
  }
  
}


void build_jcr_forest(double *x, int *nsample , int *ntree, int *maxnodes) {
  
  double bla = 4.0;
  int xdim = *nsample;
  //*ntree = bla;
  int row = 5;
  int col = 1;
  
  x[row + col*xdim] = bla;
  
  //maxnodes = 2 + ntree;
  // allocating the arrays for to store the nodes
  /*int ldaughter[ntree][maxnodes]; // 
  int rdaughter[ntree][maxnodes]; //
  int state[ntree][maxnodes]; // distinguish terminal and internal nodes
  int test_var[ntree][maxnodes]; // which variable does the node test on
  int test_thr[ntree][maxnodes]; // the value for the test
  
  build_jcr_tree(&ldaughter[2][0],maxnodes);
  */

}
