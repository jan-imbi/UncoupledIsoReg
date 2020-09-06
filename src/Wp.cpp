#include "Wp.h"



//' Minimize the entropic Wasserstein distance
//'
//'
//'
//' @param mu_vals,nu_vals numeric vector of values of discrete distribution
//' @param mu_probs,nu_probs numeric vector of probabilities of discrete distribution
//' @param p exponent of L_p distance
//'
//' @examples
//' Wp(c(1,2,3), c(1/3, 1/3, 1/3), c(2,3,4), c(1/4, 1/4, 1/2), 1)
//' @export
// [[Rcpp::export]]
double Wp(const arma::colvec& mu_vals,
              const arma::colvec& mu_probs,
              const arma::colvec& nu_vals,
              const arma::colvec& nu_probs,
              const double& p) {
  double r = 0;
  double sm = 0;
  uword iter_mu = 0;
  uword iter_nu = 0;
  uword l_mu = mu_probs.n_elem;
  uword l_nu = nu_probs.n_elem;
  colvec mu_probs_ = cumsum(mu_probs);
  colvec nu_probs_ = cumsum(nu_probs);

  for (uword i=0; i < (l_mu + l_nu - 1); ++i){
    if ((iter_mu > l_mu - 1) || (iter_nu > l_nu - 1))
      break;

    if (mu_probs_(iter_mu) < nu_probs_(iter_nu)){
      sm += (mu_probs_(iter_mu) - r) * fabs(pow(mu_vals(iter_mu) - nu_vals(iter_nu), p));
      r = mu_probs_(iter_mu);
      ++iter_mu;
    }
    else {
      sm += (nu_probs_(iter_nu) - r) * fabs(pow(mu_vals(iter_mu) - nu_vals(iter_nu), p));
      r = nu_probs_(iter_nu);
      ++iter_nu;
    }
  }
  return pow(sm, 1/p);
}
