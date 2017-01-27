#include <Rcpp.h>
using namespace Rcpp;

void reallytimeTwo(int& x){
  x = 2*x;
}

//' Multiply a number by two
//' 
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
Rcpp::List timesTwo(int x) {
   reallytimeTwo(x);
   int y = x/2;
   return Rcpp::List::create(_["double_x"]=x,
                             _["original_x"]=y);
}
