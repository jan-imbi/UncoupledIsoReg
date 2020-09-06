#ifndef __ber__
#define __ber__

// #define ARMA_NO_DEBUG
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include "RcppArmadillo.h"
#include "Rcpp.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;


arma::colvec pber(const arma::colvec& x, const double& a = -1, const double& b = 1, const double& p = 0.5);
arma::rowvec pber(const arma::rowvec& x, const double& a = -1, const double& b = 1, const double& p = 0.5);
arma::mat pber(const arma::mat& x, const double& a = -1, const double& b = 1, const double& p = 0.5);
double pber(const double& x, const double& a = -1, const double& b = 1, const double& p = 0.5);

std::vector<double> pbernoulli_custom(const arma::colvec& x, const double& a, const double& b, const double& p);
NumericVector rbernoulli_custom(const int& n, const double& a, const double& b, const double& p);

#endif // __ber__
