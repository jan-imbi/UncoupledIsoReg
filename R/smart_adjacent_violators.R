smart_isoreg <- function(x, y, D_min, D_max, V_min, V_max, tol=1e-10){
  n <- length(y)
  yf <- y
  known <- rep(F, n)
  for (i in 1:(n-1)){
    if (y[i] > y[i+1] + tol){
      yf[i] <-  y[i] + D_min
      yf[i+1] <- y[i+1] + D_max
      known[c(i, i+1)] <- c(T, T)
    }
  }

  for (i in 1:n){
    if (yf[i] + tol < V_min){
      yf[i] <- yf[i]+ D_max
      known[i] <- T
    }
    if (yf[i] >  V_max + tol ){
      yf[i] <- yf[i] + D_min
      known[i] <- T
    }
  }

  violators <- which( yf[1:(n-1)] > yf[2:n] + tol)
  last_violators <- violators
  while(length(violators) > 0){
    for (i in violators){
      if (known[i] & known[i+1]){
        stop("Something went wrong. Starting function may not be isotonic.")
      }
      if (known[i+1]){
        yf[i] <- yf[i] + D_min
        known[i] <- T
      }
      else if (known[i]){
        yf[i+1] <- yf[i+1] + D_max
        known[i+1] <- T
      }
    }
    last_violators <- violators
    violators <- which( yf[1:(n-1)] > yf[2:n] + tol)
    if (all(last_violators %in% violators)){
      break;
    }
  }
  yf
}
