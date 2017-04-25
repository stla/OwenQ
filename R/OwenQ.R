#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t number, positive or negative, possibly infinite
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (upper bound of the integral) vector of finite positive numbers, with the same length as \code{delta}
#' @param jmax,cutpoint parameters controlling the algorithm for the Owen-T function;
#' see \code{\link{OwenT}} (used only when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(nu=5, t=3, delta=2, R=100)
#' pt(q=3, df=5, ncp=2)
OwenQ1 <- function(nu, t, delta, R, jmax=50L, cutpoint=8){
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
  RcppOwenQ1(nu, t, delta, R, jmax=jmax, cutpoint=cutpoint)
}

#' @title Second Owen Q-function
#' @description Evaluates the second Owen Q-function (integral from \eqn{R} to \eqn{\infty})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t number, positive or negative, possibly infinite
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (lower bound of the integral) vector of finite positive numbers,
#' with the same length as \code{delta}
#' @param jmax,cutpoint parameters controlling the algorithm for the Owen-T function;
#' see \code{\link{OwenT}} (used only when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral
#' from \eqn{R} to \eqn{\infty}.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' # OwenQ1(nu, t, delta, R) + OwenQ2(nu, t, delta, R) = pt(t, nu, delta)
#' OwenQ1(nu=5, t=3, delta=2, R=1) + OwenQ2(nu=5, t=3, delta=2, R=1)
#' pt(q=3, df=5, ncp=2)
OwenQ2 <- function(nu, t, delta, R, jmax=50L, cutpoint=8){
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
  RcppOwenQ2(nu, t, delta, R, jmax=jmax, cutpoint=cutpoint)
}
