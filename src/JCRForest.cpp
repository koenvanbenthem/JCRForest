#include <Rcpp.h>
using namespace Rcpp;

/* all functions in this file are heavily based on the source code of the randomForest package */

void build_jcr_tree(int maxnodes) {
  
  // create vector containing which node is whose left and whose right daughter
  
  //
  
  int ncur = 0;
  
  for(int i=0; i<maxnodes-2; i++){
    
  }
}


//' Build the joint classification/regression forest
//' 
//' @param x A single integer.
//' @export return a list
// [[Rcpp::export]]
Rcpp::List build_jcr_forest(int x, int maxnodes) {
  
   int y = x/2;
   
   build_jcr_tree(maxnodes);
   
   return Rcpp::List::create(_["double_x"]=x,
                             _["original_x"]=y);
}
