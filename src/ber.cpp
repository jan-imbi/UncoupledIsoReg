#include "ber.h"

//' Discrete distribution with 2 values.
//'
//' @param n number of observations.
//' @param a first value.
//' @param b second value.
//' @param p Let X ~ pbernoulli_custom. Then p = Probability(X=a).
//' @export
// [[Rcpp::export]]
NumericVector rbernoulli_custom(const int& n, const double& a=-1, const double& b=1, const double& p=0.5){
  NumericVector ret = Rcpp::rbinom(n, 1, p);
  NumericVector::iterator it = ret.begin();
  NumericVector::iterator it_end = ret.end();
  for (; it != it_end; ++it){
    if (*it == 0.0)
      *it=a;
    else
      *it=b;
  }
  return ret;
}

//' Discrete distribution with 2 values.
//'
//' @param x observations.
//' @param a first value.
//' @param b second value.
//' @param p Let X ~ pbernoulli_custom. Then p = Probability(X=a).
//' @export
// [[Rcpp::export]]
std::vector<double> pbernoulli_custom(const arma::colvec& x, const double& a, const double& b, const double& p){
  colvec ret = colvec(x);
  colvec::iterator it = ret.begin();
  colvec::iterator it_end = ret.end();
  for (; it != it_end; ++it){
    if (*it < a)
      *it = 0;
    else if (*it >= b)
      *it = 1;
    else
      *it = p;
  }
  std::vector<double> ret_v = conv_to<std::vector<double>>::from(ret);
  return ret_v;
}



arma::colvec pber(const arma::colvec& x, const double& a, const double& b, const double& p){
  colvec ret = colvec(x);
  colvec::iterator it = ret.begin();
  colvec::iterator it_end = ret.end();
  for (; it != it_end; ++it){
    if (*it < a)
      *it = 0;
    else if (*it >= b)
      *it = 1;
    else
      *it = p;
  }
  return ret;
}

arma::rowvec pber(const arma::rowvec& x, const double& a, const double& b, const double& p){
  rowvec ret = rowvec(x);
  rowvec::iterator it = ret.begin();
  rowvec::iterator it_end = ret.end();
  for (; it != it_end; ++it){
    if (*it < a)
      *it = 0;
    else if (*it >= b)
      *it = 1;
    else
      *it = p;
  }
  return ret;
}


arma::mat pber(const arma::mat& x, const double& a, const double& b, const double& p){
  mat ret = mat(x);
  mat::iterator it = ret.begin();
  mat::iterator it_end = ret.end();
  for (; it != it_end; ++it){
    if (*it < a)
      *it = 0;
    else if (*it >= b)
      *it = 1;
    else
      *it = p;
  }
  return ret;
}


double pber(const double& x, const double& a, const double& b, const double& p){
    if (x < a)
      return 0;
    else if (x >= b)
      return 1;
    else
      return p;
}


