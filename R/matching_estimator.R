W_conjecture <- function(vec, indices){
  ret <- 0;
  for (i in 1:nrow(indices)){
    ret <- ret + abs( abs(vec[indices[i, 1]] - vec[indices[i, 2]]) - 2) / 2
  }
  ret <- ret / length(vec)
  ret
}

match_exact_Sn <- function(vec){
  n <- length(vec)
  Sn <- combinat::permn(1:n)
  best_W <- 10000;
  best_indx <- matrix(c(1:n, 1:n))
  for (sigma in Sn){
    indx <- matrix(c(1:n, sigma), ncol=2)
    W <- W_conjecture(vec, indx)
    if (W < best_W){
      best_W <- W
      best_indx <- indx;
    }
  }
  return(list(W = best_W, indices = best_indx))
}

match_exact_H <- function(vec){
  n <- length(vec)
  Sn <- combinat::permn(1:n)
  H <- list()
  for (sigma in Sn){
    if (  all( ((1:n)[sigma])[sigma] == 1:n )  ){
      H <- c(H, list(sigma));
    }
  }
  best_W <- 10000;
  best_indx <- matrix(c(1:n, 1:n))
  for (sigma in H){
    indx <- matrix(c(1:n, sigma), ncol=2)
    W <- W_conjecture(vec, indx)
    if (W < best_W){
      best_W <- W
      best_indx <- indx;
    }
  }
  return(list(W = best_W, indices = best_indx))
}


match_sequentially_H <- function(vec){
  ret <- matrix(rep(0, times=2*length(vec)), nrow=length(vec));
  indices <- 1:length(vec);
  mat_size <- 0;
  while (length(indices)>0){
    idx <- indices[1];
    mat_size <- mat_size+1;

    best_diff <- 1;
    best_idx <- idx;
    for (i in indices){
      diff <- abs(abs(vec[idx] - vec[i]) -2) / 2;
      if (diff < best_diff){
        best_idx <- i;
        best_diff <- diff;
      }
    }
    ret[mat_size, ] <- c(idx, best_idx)
    indices <- indices[ !(indices %in% c(idx, best_idx)) ]
  }
  ret <- ret[1:mat_size,]

  real_ret <- matrix(rep(0, 2*length(vec)), ncol=2);
  iter <- 0;
  for (i in 1:nrow(ret)){
    iter <- iter+1
    if (ret[i, 1] == ret[i, 2]){
      real_ret[iter, ] <- ret[i, ]
    }
    else{
      real_ret[iter, ] <- ret[i, ]
      iter <- iter+1
      real_ret[iter, ] <- rev(ret[i, ])
    }
  }
  real_ret <- real_ret[order(real_ret[,1]),]
  return(real_ret)
}

match_sequentially_Sn <- function(vec){
  ret <- matrix(rep(0, times=2*length(vec)), nrow=length(vec));
  indices_one <- 1:length(vec);
  indices_two <- 1:length(vec);
  for (idx_one in indices_one){
    best_diff <- 10000;
    best_idx <- indices_two[1];
    for (i in indices_two){
      diff <- abs(abs(vec[idx_one] - vec[i])- 2) / 2;
      if (diff < best_diff){
        best_idx <- i;
        best_diff <- diff;
      }
    }
    indices_two <- indices_two[indices_two != best_idx]
    ret[idx_one, ] <- c(idx_one, best_idx)
  }
  return(ret)
}


match_random_Sn <- function(vec){
  ret <- matrix(rep(0, times=2*length(vec)), nrow=length(vec));
  indices_one <- 1:length(vec);
  indices_one <- sample(indices_one, size = length(vec), replace = F)
  indices_two <- 1:length(vec);
  indices_two <- sample(indices_two, size = length(vec), replace = F)
  for (idx_one in indices_one){
    best_diff <- 10000;
    best_idx <- indices_two[1];
    for (i in indices_two){
      diff <- abs(abs(vec[idx_one] - vec[i])- 2) / 2;
      if (diff < best_diff){
        best_idx <- i;
        best_diff <- diff;
      }
    }
    indices_two <- indices_two[indices_two != best_idx]
    ret[idx_one, ] <- c(idx_one, best_idx)
  }
  return(ret)
}


minimizing_mu_match <- function(vec, indices){
  ret <- numeric(length(vec));
  for (i in 1:nrow(indices)){
    ret[i] <- (vec[indices[i, 1]] +  vec[indices[i, 2]])/2.0
  }
  ret
}
