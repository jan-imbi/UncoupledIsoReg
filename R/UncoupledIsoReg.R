#' @details
#' \code{UncoupledIsoReg} is a package for uncoupled isotonic regression.
#' @keywords internal
"_PACKAGE"

# .onUnload <- function (libpath) {
#   library.dynam.unload("UncoupledIsotonicRegression.dll", libpath)
#   unloadNamespace("UncoupledIsotonicRegression")
#   }

.onUnload <- function (libpath) {
  library.dynam.unload("UncoupledIsoReg", libpath)
}

#' @useDynLib UncoupledIsoReg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
