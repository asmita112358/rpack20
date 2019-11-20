// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]



double d_cov(  arma::vec x, arma::vec y, double N)
  {
  int n = x.size();
  int Upper=100;
  float c=pow(3.1415,2)/2;
  
  arma::colvec  u=linspace(1,100,100);
  
  arma::colvec v = 1/c*1/pow(u-0.5,2);
  
  arma::colvec p = v/sum(v);
  arma::colvec k = Rcpp::RcppArmadillo::sample(u,N,TRUE,p);
  arma::colvec kv= (k-0.5)*3.1415;
  arma::mat U= x*trans(kv);
  arma::mat V= y*trans(kv);
  arma::mat phix= sqrt(2/N)*sin(U);
  arma::mat phiy= sqrt(2/N)*sin(V);
  arma::mat Kx= phix*trans(phix);
  arma::mat Ky= phiy*trans(phiy);
  Kx.diag(0).fill(0);
  Ky.diag(0).fill(0);
  arma::mat rx=sum(Kx);
  arma::mat ry=sum(Ky);
  double t3=2*as_scalar(rx*trans(ry))/(n-2);
  double t1=as_scalar(trace(Kx*Ky));
  double t2=accu(rx)*accu(ry)/((n-1)*(n-2));
  double d=(t1+t2-t3)/((n-3)*n);
  return d;
  //return List::create( Named("Beta") =t3);
  
}