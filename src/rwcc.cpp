#include <Rcpp.h>
using namespace Rcpp;

//' @title A Metropolis sampler using Rcpp 
//' @description A Metropolis sampler using Rcpp to generate random variables from a standard Cauchy distribution
//' @param sigma the parameter of proposal distribution Normal
//' @param x0 the initial value of sample
//' @param N the number of samples
//' @return a random sample of size N and the number of reject 
//' @examples
//' \dontrun{
//' rwC <- rwcc(0.5,25,2000)
//' print(rwC$x)
//' print(rwC$k)
//' }
//' @export
// [[Rcpp::export]]
List rwcc(double sigma, double x0, int N) {
  Rcpp::NumericVector x(N) ;
  x[0] = x0;
  Rcpp::NumericVector u ;
  u= runif(N);
  int k = 0;
  for (int i=1;i<N; i+=1) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1]))))
    {x[i] = y[0] ;
    }else  {
      x[i] = x[i-1];
      k = k + 1;
    } }
  return List::create(Named("x") =x, Named("k") = k);
}
