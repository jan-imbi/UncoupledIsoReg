#' Convert a discrete probability distribution to an isotonic function
#'
#'
#' @param mu A named list containing two vectors of the same length, \code{`vals`} and \code{`probs`}.
#' mu represents a discrete probability distribution.
#' @param n integer number of equidistant points on [min(mu$vals), max(mu$vals)] which will be the support of the
#' returned isotonic function
#' @param tol tolerance for shifting cumulative probabilities
#' @return Vector representing isotonic function values
#' @export
#' @importFrom utils head tail
measure_to_iso <- function(mu, n=NULL, tol = .Machine$double.eps * 1e4) {
  if(!all.equal(sum(mu$probs), 1)){
    stop("Probabilities of mu dont add up to one.")
  }
  if (is.null(n)){
    n <- 3*length(mu[["probs"]])
  }

  ret <- double(n)
  e_cdf <- cumsum(mu$probs)
  tolerance_fix_n <- floor(length(e_cdf)/2)
  e_cdf[1:tolerance_fix_n] <- e_cdf[1:tolerance_fix_n] - tol
  e_cdf[tolerance_fix_n:length(e_cdf)] <-  e_cdf[tolerance_fix_n:length(e_cdf)] + tol
  for (i in seq.int(length.out = n)){
    ret[i] <- mu$vals[which(e_cdf >= i/n)[1]]
  }
  return(ret)
}


#' Convert a discrete probability distribution to a smooth mostly-isotonic function
#'
#' @param mu A named list containing two vectors of the same length, \code{`vals`} and \code{`probs`}.
#' mu represents a discrete probability distribution.
#' @param n integer number of equidistant points on [min(mu$vals), max(mu$vals)] which will be the support of the
#' returned isotonic function
#' @param tol tolerance for shifting cumulative probabilities
#' @param kernel the kernel to be used. Can be abbreviated. Type \code{?ksmooth} for details.
#' @param bandwidth the bandwidth. The kernels are scaled so that their
#' quartiles (viewed as probability densities) are at +/- 0.25*bandwidth. Type \code{?ksmooth} for details.
#' @return Vector representing isotonic function values
#' @export
#' @importFrom stats ksmooth
#' @importFrom utils head tail
measure_to_smooth_iso <- function(mu, n, kernel="normal", bandwidth=.1, tol = .Machine$double.eps * 10){
  if(!all.equal(sum(mu$probs), 1)){
    stop("Probabilities of mu dont add up to one.")
  }
  if (is.null(n)){
    n <- 3*length(mu[["probs"]])
  }

  ret <- double(n)
  e_cdf <- cumsum(mu$probs)
  tolerance_fix_n <- floor(length(e_cdf)/2)
  e_cdf[1:tolerance_fix_n] <- e_cdf[1:tolerance_fix_n] - tol
  e_cdf[tolerance_fix_n:length(e_cdf)] <-  e_cdf[tolerance_fix_n:length(e_cdf)] + tol
  for (i in seq.int(length.out = n)){
    ret[i] <- mu$vals[which(e_cdf >= i/n)[1]]
  }
  X <- seq(from=head(mu$vals, 1), to=tail(mu$vals, 1), length.out =n)
  ksmooth(X, ret, kernel=kernel, bandwidth=bandwidth, x.points=X)$y
}



project <- function(from, to){
  ret <- double(length(from))
  for (i in seq_along(from)){
    ind <- which.min(abs(from[[i]] - to))
    ret[i] <- to[ind]
  }
  ret
}

pushforward <- function(mu, f, ...){
  ret <- list()
  vals <- f(mu$vals, ...)
  ret$vals <- sort(unique(vals))
  ret$probs <- numeric(length(ret$vals))
  for (i in seq_along(ret$prob)){
    ret$probs[[i]] <- sum(vals==ret$vals[[i]]);
  }
  ret$probs <- ret$probs / sum(ret$probs)
  ret
}

















