#ifndef __Wp__
#define __Wp__

#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

double Wp(const colvec& mu_vals, const colvec& mu_probs_, const colvec& nu_vals, const colvec& nu_probs_, const double& p);

#endif // __Wp__
