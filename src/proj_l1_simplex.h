#ifndef __proj_l1_simplex__
#define __proj_l1_simplex__

#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"
using namespace arma;
using namespace Rcpp;

arma::colvec laurentCondat(const arma::colvec& y);

#endif // __proj_l1_simplex__
