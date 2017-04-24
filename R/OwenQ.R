#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t finite number, positive or negative
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (upper bound of the integral) vector of finite positive numbers, with the same length as \code{delta}
#' @param jmax,cut.point passed to \code{\link{OwenT}} (when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(nu=4, t=3, delta=2, R=100)
#' ptOwen(q=3, nu=4, delta=2)
OwenQ1 <- function(nu, t, delta, R, jmax=50L, cut.point=8){
  if(length(delta) != length(R)){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0)){
    stop("`R` must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(R))){
    stop("`R` must be finite.")
  }
  if(any(is.infinite(delta))){
    stop("`delta` must be finite.")
  }
  RcppOwenQ1(nu, t, delta, R, jmax=jmax, cutpoint=cut.point)
}
