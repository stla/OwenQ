#' @title Owen's equality 11
#' @description Evaluates the Owen cumulative distribution function in the 4th case.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t1,t2 two finite numbers, positive or negative
#' @param delta1,delta2 two vectors of finite numbers, with the same length
#' @param jmax,cutpoint parameters controlling the algorithm for the Owen-T function;
#' see \code{\link{OwenT}} (used only when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' # Wolfram integration gives 0.018689824158
#' pOwen4(nu=5, t1=2, t2=1, delta1=3, delta2=2)
pOwen4 <- function(nu, t1, t2, delta1, delta2, jmax=50L, cutpoint=8){
  J <- length(delta1)
  if(J != length(delta1)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2)){
    stop("`delta1` must be >`delta2`.")
  }
  if(t1<t2){
    stop("`t1` must be >`t2`.")
  }
  if(is.infinite(t1) || is.infinite(t2)){
    stop("`t1` and `t2` must be finite.")
  }
  if(nu == Inf){
    return(pmax(0, pnorm(t2, mean=delta2)-pnorm(t1, mean=delta1)))
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(t1==t2){
    return(ptOwen(t2, nu, delta2, jmax, cutpoint)-ptOwen(t2, nu, delta1, jmax, cutpoint))
  }
  RcppOwenCDF4(nu, t1, t2, delta1, delta2, jmax, cutpoint)
}
