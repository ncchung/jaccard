#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double kahanSum(NumericVector x) {
  double sum = 0.0;
  double y,t;
  double c = 0.0;
  NumericVector::iterator it;
  
  for( it = x.begin(); it != x.end(); ++it ){
    y = *it - c;
    t = sum + y;
    c = (t-sum)-y;
    sum = t;
  }
  return sum;
}


// [[Rcpp::export]]
NumericVector kahanCumSum(NumericVector x) {
  double 	sum = 0.0;
  double 	y,t;
  double 	c = 0.0;
  int 		N = x.size();
  NumericVector res(N);

  for( int i=0; i<N; ++i){
    y = x[i] - c;
    t = sum + y;
    c = (t-sum)-y;
    sum = t;
    res[i] = sum;
  }
  return res;
}


