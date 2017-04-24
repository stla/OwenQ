#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
OwenT01 <- function(h, a, jmax=50L, cutpoint=8) {
  if(isNotPositiveInteger(jmax)){
    stop("`jmax` must be an integer >=1.")
  }
  if(cutpoint <= 0){
    stop("`cutpoint` must be a strictly positive number.")
  }
  if(a<0 || a>1){
    stop("`a` must be a number between 0 and 1.")
  }
  if(h<0){
    stop("`h` must be positive.")
  }
  RcppOwenT01(h, a, jmax, cutpoint)
}

#' @title Owen T-function
#' @description Evaluates the Owen T-function.
#' @param h numeric scalar
#' @param a numeric scalar
#' @param jmax integer controlling the number of terms of the
#' series expansion; see Details
#' @param cutpoint positive number, the cut point in the algorithm;
#' see Details
#' @return A number between 0 and 1.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @export
#' @details If \eqn{0\lea\le1}, and \eqn{0\leh\lec}, where \eqn{c} is
#' the cut point, a series expansion is used.
#' It is truncated after the \code{jmax}-th term.
#' If \eqn{0\lea\le1}, and \eqn{h>c}, an asymptotic approximation is used.
#' Otherwise, the properties of the Owen T-function are exploited
#' to come down to the case \eqn{0\lea\le1}.
#' See the reference for more information.
#' @references
#' Owen, D. B. (1956).
#' Tables for computing bivariate normal probabilities.
#' \emph{Ann. Math. Statist.} \bold{27}, 1075-1090.
#' @examples
#' integrate(function(x) pnorm(1+2*x)^2*dnorm(x), lower=-Inf, upper=Inf)
#' pnorm(1/sqrt(5)) - 2*OwenT(1/sqrt(5), 1/3)
OwenT <- function (h, a, jmax = 50L, cutpoint = 8)
{
  if(isNotPositiveInteger(jmax)){
    stop("`jmax` must be an integer >=1.")
  }
  if(cutpoint <= 0){
    stop("`cutpoint` must be a strictly positive number")
  }
  if (!is.vector(a) || length(a) > 1L)
    stop("'a' must be a vector of length 1")
  if (!is.vector(h))
    stop("'h' must be a vector")
  if (is.na(a))
    stop("parameter 'a' is NA")
  if (is.infinite(h)){
    return(0)
  }
  return(RcppOwenT(h, a, jmax, cutpoint))
}
