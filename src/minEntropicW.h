#ifndef __minEntropicW__
#define __minEntropicW__

#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

Rcpp::List minimize_entropic_W(
    const arma::colvec& Y, const double& V_min, const double& V_max,
    const arma::uword& N,
    Rcpp::Nullable<arma::colvec> A,
    Rcpp::Nullable<arma::colvec> AV,
    Rcpp::Nullable<arma::mat> P_D,
    Rcpp::Nullable<arma::colvec> muStart,
    const arma::uword& maxIter, const arma::uword& sinkhornIter,
    const double& eps,
    const double& gammaStart, const double& p,
    const double& sinkhornTol,
    const double& gradDescTol,
    const std::string& noise,
    const double& sigma,
    const double& D_min, const double& D_max, const double& p_ber,
    const bool& fastSinkhorn,
    const bool& pushforwardStart,
    const bool& suppressOutput,
    const double& WThreshold);

Rcpp::List min_entropic_W_calc(
      const arma::colvec& Y, const arma::colvec& A, const arma::colvec& AV,
      const arma::mat& P,
      const arma::colvec& muStart,
      const arma::uword& maxIter,
      const arma::uword& sinkhornIter,
      const double& eps,
      const double& gammaStart,
      const double& p,
      const double& sinkhornTol,
      const double& gradDescTol,
      const bool& fastSinkhorn,
      const bool& suppressOutput,
      const double& WThreshold);

#endif // __minEntropicW__
