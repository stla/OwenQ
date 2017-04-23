#' @title Owen's equality 11
#' @description Evaluates the Owen cumulative distribution function in the 4th case.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t1,t2 two finite numbers, positive or negative
#' @param delta1,delta2 two vectors of finite numbers, with the same length
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
pOwen4 <- function(nu, t1, t2, delta1, delta2){
  J <- length(delta1)
  if(J != length(delta1)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2)){
    stop("`delta1` must be >`delta2`.")
  }
  if(any(t1<=t2)){
    stop("`t1` must be >`t2`.")
  }
  if(any(is.infinite(t1) | is.infinite(t2))){
    stop("`t1` and `t2` must be finite.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  RcppOwenQ1(nu, t2, delta2, R, jmax=8L, cutpoint=50) -
    RcppOwenQ1(nu, t1, delta1, R, jmax=8L, cutpoint=50)
}
