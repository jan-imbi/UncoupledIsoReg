#include "minEntropicW.h"
#include "sinkhorn.h"
#include "Wp.h"
#include "ber.h"
#include "proj_l1_simplex.h"

arma::colvec calculateNearestNeighbourDistribution_cpp(arma::colvec outputSpace_vals, const arma::colvec& comparativeSpace_vals){
  colvec ret = colvec(outputSpace_vals.n_elem, fill::zeros);
  colvec m = colvec(outputSpace_vals.n_elem);
  uword indx =0;
  for (colvec::const_iterator it = comparativeSpace_vals.begin(); it != comparativeSpace_vals.end(); ++it){
    m = arma::abs(outputSpace_vals - *it);
    indx = m.index_min();
    ret[indx] += 1;
  }
  // ret = laurentCondat(ret);
  ret = ret/accu(ret);
  return ret;
}

Rcpp::List min_entropic_W_calc(
    const arma::colvec& Y,
    const arma::colvec& A,
    const arma::colvec& AV,
    const arma::mat& P,
    const arma::colvec& muStart,
    const arma::uword& maxIter = 100,
    const arma::uword& minIter = 50,
    const arma::uword& sinkhornIter = 100,
    const double& eps = 0.01,
    const double& gammaStart = 0.01,
    const double& p=1,
    const double& sinkhornTol = 1e-12,
    const double& gradDescTol = 1e-12,
    const bool& fastSinkhorn = true,
    const bool& suppressOutput = false,
    const double& WThreshold = 0){

  const colvec& Y_sorted = sort(Y);
  const arma::colvec& a=A;
  const arma::colvec& av=AV;
  colvec mu_hat = muStart;

  mat C = mat(a.n_elem, Y_sorted.n_elem);
  for (uword i=0; i<a.n_elem; ++i){
    for (uword j=0; j<Y_sorted.n_elem; ++j){
      C(i,j) = fabs(pow(a(i)-Y_sorted(j), p));
    }
  }

  colvec pi_hat = colvec(Y_sorted.n_elem);
  pi_hat.fill(1.0/Y_sorted.n_elem);

  colvec nu_probs = P * mu_hat;

  double L_new;
  double W_new;
  field<colvec> lsnk;

  const mat P_t = P.t();
  colvec grad = colvec(mu_hat.n_elem);
  colvec grad_last = colvec(mu_hat.n_elem);
  colvec mu_last = colvec(mu_hat.n_elem);
  uword iter=0;

  arma::colvec best_mu=mu_hat;
  arma::colvec mu_threshold;
  double best_W = std::numeric_limits<double>::max();
  double best_L = std::numeric_limits<double>::max();
  double gamma = gammaStart;

  std::vector<double> L_vals;
  std::vector<double> W_vals;
  L_vals.reserve(maxIter);
  W_vals.reserve(maxIter);

  bool fast_fail=false;


  for (iter=0; iter<maxIter; ++iter){
    if (fastSinkhorn){
      lsnk = Subgradient(nu_probs, pi_hat, C, eps, sinkhornIter, sinkhornTol);
      if (!(lsnk(0).is_finite()) | (norm(lsnk(0)) < 1e-10)){
        fast_fail=true;

        if (!suppressOutput ){
          Rcout << "Fast Sinkhorn failed in iteration "<< iter <<". Entering stabilized fallback." << std::endl;
          Rcout << norm(lsnk(0)) << std::endl;
        }
      }
    }
    if (!fastSinkhorn || fast_fail){
      if (fast_fail){
        lsnk = logStabSinkhorn_cpp(nu_probs, pi_hat, C, eps, sinkhornIter/3, sinkhornTol);
        fast_fail=false;
        if (!suppressOutput ){
          Rcout << "Finished fallback mode." << std::endl;
        }

      }
      else{
        lsnk = logStabSinkhorn_cpp(nu_probs, pi_hat, C, eps, sinkhornIter, sinkhornTol);
      }
    }

    grad = P_t * lsnk(0);
    if (iter > 0){
      gamma = fabs(dot((mu_hat - mu_last), (grad - grad_last)) / pow( norm(grad - grad_last), 2));
    }

    mu_last = mu_hat;
    grad_last = grad;

    mu_hat = mu_last - gamma*grad;
    mu_hat = laurentCondat(mu_hat);
    nu_probs = P * mu_hat;

    double L_last = L_new;
    L_new=as_scalar(lsnk(1));
    L_vals.push_back(L_new);
    W_new = Wp(a, P*mu_last, Y_sorted, pi_hat, p);
    W_vals.push_back(W_new);

    if (W_new < best_W){
      best_W = W_new;
      best_L = L_new;
      best_mu = mu_last;
    }

    if (W_new < WThreshold){
      mu_threshold = mu_last;
      break;
    }

    if ((fabs(L_new - L_last) < gradDescTol) &&(iter > minIter)){
      break;
    }
  }

  if (!suppressOutput ){
    if (iter == maxIter){
      Rcout << "Maximum number iterations reached (" << iter << ")." << std::endl;
    }
    else{
      Rcout << "Convergence criterion met after " << iter << " iterations." << std::endl;
    }
  }
  L_vals.shrink_to_fit();
  W_vals.shrink_to_fit();

  return Rcpp::List::create(Rcpp::Named("vals") = conv_to<std::vector<double>>::from(av),
                            Rcpp::Named("probs") = conv_to<std::vector<double>>::from(best_mu),
                            Rcpp::Named("A") = conv_to<std::vector<double>>::from(a),
                            Rcpp::Named("P") = P,
                            Rcpp::Named("Best_W_val") = best_W,
                            Rcpp::Named("Best_L_val") = best_L,
                            Rcpp::Named("L_vals") = L_vals,
                            Rcpp::Named("W_vals") = W_vals,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("last_mu") = conv_to<std::vector<double>>::from(mu_hat),
                            Rcpp::Named("mu_threshold") = conv_to<std::vector<double>>::from(mu_threshold)
  );
}


//' Minimize the entropic Wasserstein distance
//'
//'
//'
//' @param Y numeric vector of observations
//' @param A numeric vector of domain values where you expect most of the Y values to lie in
//' @param AV numeric vector of domain values which you allow the functions in the space you minimize over to take
//' @param P_D Jacobian of mu -> W_2(mu*D, pi_hat)
//' @param muStart starting values for minimization
//' @param maxIter maximum iteration number for gradient descent algorithm
//' @param minIter minimum iteration number for gradient descent algorithm
//' @param sinkhornIter maximum Sinkhorn iterations
//' @param eps entropic regularization constant
//' @param gammaStart starting gamma value for gradient descent algorithm. Only used in first descent step.
//' @param p exponent of L_p norm
//' @param sinkhornTol tolerance for stopping criterion in Sinkhorn algorithm
//' @param gradDescTol tolerance for stopping criterion in gradient descent algorithm based on euclidian distance to last
//' vector
//' @param fastSinkhorn logical, controls whether to use the fast Sinkhorn algorithm
//' @param pushforwardStart if TRUE, starts with nearest neighbour distribution to pi_hat instead of rep(1/n)
//' @param suppressOutput suppress output messages?
//' @param WThreshold tolerance for stopping criterion for gradient descent algorithm based on the Wasserstein distance
//'
//' @examples
//' library(tidyverse)
//' n <- 1000
//' x <- seq(0, 1, length.out = n)
//' m <- function(x) (2*(x- 0.5))^3
//' Y_no_error <-  m(x)
//' varepsilon <- rbernoulli_custom(n, a= -0.3, b= 0.3, p=0.5)
//' Y <- (Y_no_error + varepsilon) %>% sample(n)
//' dat <- tibble(x=x, Y=Y, Y_no_error = Y_no_error)
//' N <- round(sqrt(n))
//' A <- seq(-1.3, 1.3, length.out = N)
//' stepsize <- (A[2]-A[1]) / 2
//' A_V <-  A[which(-1 - stepsize <= A & A <= 1 + stepsize)]
//' p_ber <- function(x) pbernoulli_custom(x, a = -0.3, b = 0.3, p=0.5)
//' P <- matrix(rep(0, times = (N) * length(A_V)), nrow = N)
//' P[1, ] <- p_ber(A[2] - stepsize - A_V)
//' for (i in 2:(N-1)) {
//'     P[i,] <- p_ber(A[i + 1] - stepsize - A_V) - p_ber(A[i] - stepsize - A_V)
//' }
//' P[N, ] <- 1 - p_ber(A[N] - stepsize - A_V)
//'
//' l <-
//'    minimize_entropic_W(Y = Y,
//'                        A = A,
//'                        AV = A_V,
//'                        P_D = P,
//'                        suppressOutput = FALSE)
//' mu_hat <- list()
//' mu_hat[["vals"]] <- l$vals
//' mu_hat[["probs"]] <- l$probs
//' g_hat <- measure_to_smooth_iso(mu_hat, n)
//' dat <- bind_cols(dat, g_hat = g_hat)
//' ggplot(dat, aes(x=x, y=Y_no_error)) +
//'   geom_line(col="red") +
//'   scale_y_continuous("") +
//'   geom_line(aes(y=g_hat), col="blue")
//'
//' @import tidyverse
//' @export
// [[Rcpp::export]]
Rcpp::List minimize_entropic_W(
    const arma::colvec& Y,
    Rcpp::Nullable<arma::colvec> A = R_NilValue,
    Rcpp::Nullable<arma::colvec> AV = R_NilValue,
    Rcpp::Nullable<arma::mat> P_D = R_NilValue,
    Rcpp::Nullable<arma::colvec> muStart = R_NilValue,
    const arma::uword& maxIter = 100,
    const arma::uword& minIter = 50,
    const arma::uword& sinkhornIter = 100,
    const double& eps = 0.01,
    const double& gammaStart = .05,
    const double& p=1,
    const double& sinkhornTol = 1e-12,
    const double& gradDescTol = 1e-12,
    const bool& fastSinkhorn = true,
    const bool& pushforwardStart=false,
    const bool& suppressOutput = false,
    const double& WThreshold = 0){

  colvec a;
  colvec av;
  colvec mu_hat;
  mat P;

  a = Rcpp::as<arma::colvec>(A);
  av = Rcpp::as<arma::colvec>(AV);
  P = Rcpp::as<arma::mat>(P_D);

  if (muStart.isNull()){
    colvec start_;
    if (pushforwardStart){
      start_= calculateNearestNeighbourDistribution_cpp(av, sort(Y));
    }
    else{
      start_ = colvec(av.n_elem, fill::ones);
      start_ = start_ / (double)start_.n_elem;
    }

    return min_entropic_W_calc(Y, a, av, P, start_, maxIter, minIter, sinkhornIter, eps,
                               gammaStart, p, sinkhornTol, gradDescTol, fastSinkhorn,
                               suppressOutput, WThreshold);
  }
  else{
    return min_entropic_W_calc(Y, a, av, P, as<arma::colvec>(muStart),
                               maxIter, minIter, sinkhornIter, eps, gammaStart, p, sinkhornTol,
                               gradDescTol, fastSinkhorn,
                               suppressOutput, WThreshold);
  }
}

