#include "proj_l1_simplex.h"

// [[Rcpp::export]]
arma::colvec laurentCondat(const arma::colvec& y) {
  if ((sum(y)==1)&&(all(y>=0))){
    return(y);
  }
  int length = y.n_elem;
  arma::colvec v = arma::colvec(length);
  const double& alpha=1;
  int vend=1;
  int vendhelp;
  int vstart=-1;
  int vhelp=0;
  double rho = (v[0] = y[0]) - alpha;
  int i=1;
  for (; i<length; i++){
    if (y[i]>rho) {
      if ((rho+= ( (v[vend]=y[i]) -rho)/(vend-vstart)) <= y[i]-alpha) {
        rho=y[i]-alpha;
        vstart=vend-1;
      }
      vend++;
    }
  }
  if (vstart>=0) {
    for(i=vstart; i>=0; i--)
      if (v[i]>rho){
        double pls = ((v[vstart--]=v[i])-rho);
        pls = pls/ (vend-vstart-1);
        rho+=pls;
      }
  }
  vhelp=-2;
  while (vhelp<i) {
    vendhelp = vend;
    for (i=vhelp=vstart+1; i<vendhelp; i++)
      if (v[i]>rho)
        v[vhelp++]=v[i];
      else
        rho+=(rho-v[i])/((--vend)-(vstart+1));
  }
  for (i=0; i<length; i++){
    v[i]=(y[i]>rho ? y[i]-rho : 0.0);
  }
  return(v);
}
