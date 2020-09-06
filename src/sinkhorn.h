#ifndef __sinkhorn__
#define __sinkhorn__

#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

template <typename T>
double softMinStabilized_cpp(const T& z, const double& eps);

arma::field<arma::colvec> logStabSinkhorn_cpp(const arma::colvec& a_, const arma::colvec& b,const arma::mat& costMat_, const double& eps = 1e-2,
                                              const arma::uword& maxIter=200, const double& tol=1e-3);

arma::field<arma::colvec> Subgradient(const arma::colvec& a_, const arma::colvec& b, const arma::mat& M,
                                      const double& eps, const arma::uword& maxIter,
                                      const double& tolerance);

#endif // __sinkhorn__
