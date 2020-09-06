#include "sinkhorn.h"

template <typename T>
double softMinStabilized_cpp(const T& z, const double& eps){
  const double& zmin = z.min();
  return zmin - eps * log(accu(exp(- (z - zmin)/eps)));
}

// [[Rcpp::export]]
arma::field<arma::colvec> logStabSinkhorn_cpp(const arma::colvec& a_, const arma::colvec& b,const arma::mat& costMat_, const double& eps,
                          const arma::uword& maxIter, const double& tol) {
  bool ZeroValues = false;
  arma::uvec ids;
  arma::uword n_a_orig = a_.n_elem;
  if(all(a_)==false){
    ZeroValues = true;
    ids = find(a_ > 0); // Find indices with support
  }

  const arma::mat& costMat = (ZeroValues) ? costMat_.rows(ids) : costMat_;
  const arma::colvec& a = (ZeroValues) ? a_.rows(ids) : a_;

  const arma::colvec& eps_loga = eps * log(a);
  const arma::rowvec& eps_logb = conv_to<arma::rowvec>::from(eps * log(b));
  const arma::rowvec& eps_logb_colvec = (eps * log(b));
  arma::colvec f_new = arma::colvec(a.n_elem, fill::zeros);
  arma::colvec f_old = arma::colvec(a.n_elem);
  arma::rowvec g_new = arma::rowvec(b.n_elem, fill::zeros);
  arma::rowvec g_old = arma::rowvec(b.n_elem);

  arma::colvec ones_b = arma::colvec(a.n_elem, fill::ones);
  arma::colvec ones_a = arma::colvec(b.n_elem, fill::ones);

  const arma::uword ncolC = costMat.n_cols;
  const arma::uword nrowC = costMat.n_rows;
  mat S = arma::mat(size(costMat));

  for (uword iter = 0; iter<maxIter; ++iter){
    f_old = f_new;
    g_old = g_new;
    S = costMat - repmat(mat(f_old), 1, ncolC) - repmat(mat(g_old), nrowC, 1);


    // Rcpp::Rcout << (exp(-S / eps)*ones_a - a)<< endl;
    // if  ((norm(exp(-S / eps)*ones_a - a)  < tol)&&(norm(exp(-S.t() / eps)*ones_b - b)  < tol)){
    // if  ((norm(exp(-S.t() / eps)*ones_b - b)  < tol)){
    if  ((norm(exp(-S / eps)*ones_a - a)  < tol)){
      break;
    }
    for (uword i=0; i<S.n_rows; ++i){
      f_new(i) = softMinStabilized_cpp(S.row(i), eps);
    }
    f_new = f_new + f_old + eps_loga;

    S = costMat - repmat(mat(f_new), 1, ncolC) - repmat(mat(g_old), nrowC, 1);

    for (uword i=0; i<S.n_cols; ++i){
      g_new(i) = softMinStabilized_cpp(S.col(i), eps);
    }
    g_new = g_new + g_old + eps_logb;
  }

  arma::field<arma::colvec> ret(2);
  if (ZeroValues){
    arma::colvec ret_zero = arma::colvec(n_a_orig, arma::fill::zeros);
    ret_zero(ids) = f_new - accu(f_new)/ (double)f_new.n_elem;
    ret(0) = ret_zero;
    // ret(1) = conv_to<arma::colvec>::from(g_new);
  }
  else{
    ret(0) = f_new;
    // ret(1) = conv_to<arma::colvec>::from(g_new);
  }
  ret(1)= -eps*accu(exp(-S / eps)) + accu(f_old % a) + accu(g_old % b);
  return ret;
}


// This is a modified version of the Subgradient method from the Barycenter package, see
// https://cran.r-project.org/web/packages/Barycenter/index.html
// [[Rcpp::export]]
arma::field<arma::colvec> Subgradient(const arma::colvec& a_, const arma::colvec& b, const arma::mat& M,
                                      const double& eps, const arma::uword& maxIter = 100,
                                      const double& tolerance = 1e-14) {

  // Transforming Input, i.e. calculating the kernel

  const double lambda = 1/eps;

  arma::mat K = exp(-lambda*M);

  bool ZeroValues = false;
  arma::uvec ids;
  arma::mat Ktemp;
  arma::mat Mtmp = M;

  uword a_n = a_.n_elem;

  arma::colvec a;
  if(all(a_)==false){
    ZeroValues = true;
    Ktemp = K;  //save original matrix K to compute the transportplan at the end
    ids = find(a_ > 0); // Find indices with support
    K = K.rows(ids);
    a = a_.rows(ids);
    Mtmp = M.rows(ids);
  }
  else{
    a=a_;
  }
  arma::mat ainvK = diagmat(1/a) * K;


  //Initialize u, v and the vector next for the stopping criterion
  arma::vec u(a.n_rows);
  u = u.ones()/a.n_rows;
  arma::vec next;

  arma::uword i=0;
  for(; i<maxIter;i++){


    u = 1 / (ainvK * (b / (K.t()*u)));  //Sinkhorn`s update

    // The stopping criterion is checked every 20th step, i.e. if u and the next update of u (called next) differ in 2-norm only by a tolerance
    if(i % 20 == 0){
      next = 1 / (ainvK * (b / (K.t()*u)));
      double Criterion = norm(abs(next-u));
      if ((Criterion < tolerance)){
        u = next;
        break;
      }
      else{
        u = next;
      }
    }
  }

  arma::colvec v = b / (K.t()*u);
  arma:: mat P = arma::diagmat(u) * K * arma::diagmat(v);

  double L_val = accu(P % Mtmp) + eps * accu(P % (P - 1));

  arma::field<arma::colvec> ret_field(2);

  if(ZeroValues){
    arma::colvec alpha = (1 / lambda) * log(u) - accu(log(u))/(lambda*(u.n_rows));
    // arma::colvec alpha = (1 / lambda) * log(u);
    arma::colvec ret = arma::colvec(a_n, fill::zeros);
    ret(ids) = alpha;
    // ret =  ret - accu(ret)/ (double)ret.n_elem;
    ret_field(0)=ret;
  }
  else{
    arma::colvec alpha = (1 / lambda) * log(u) - accu(log(u))/(lambda*(u.n_rows));
    ret_field(0)=alpha;
  }
  ret_field(1)=L_val;
  return(ret_field);
}
















