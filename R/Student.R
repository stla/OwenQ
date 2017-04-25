#' @title Student CDF with integer number of degrees of freedom
#' @description Cumulative distribution function of the noncentrel Student
#' distribution with an integer number of degrees of freedom.
#' @param q quantile
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param delta numeric vector of noncentrality parameters
#' @param jmax,cutpoint parameters controlling the algorithm for the Owen-T function;
#' see \code{\link{OwenT}} (used only when \code{nu} is odd)
#' @return Numeric vector, the CDF evaluated at \code{q}.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of
#' freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' ptOwen(2, 3) - pt(2, 3)
#' ptOwen(2, 3, delta=1) - pt(2, 3, ncp=1)
ptOwen <- function(q, nu, delta=0, jmax=50L, cutpoint=8){
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(q) || is.infinite(delta)){
    stop("Parameters must be finite.")
  }
  RcppOwenStudent(q, nu, delta, jmax, cutpoint)
}

