#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
OwenT01 <- function(h, a, jmax=50L, cut.point=8) {
  if(isNotPositiveInteger(jmax)){
    stop("`jmax` must be an integer >=1.")
  }
  if(cut.point <= 0){
    stop("`cut.point` must be a strictly positive number")
  }
  if(a<0 || a>1){ # ok pour a = 0 ?
    stop("`a` must be a number between 0 and 1")
  }
  RcppOwenT01(h, a, jmax, cut.point)
}

#' @title Owen T-function
#' @description Evaluates the Owen T-function.
#' @param h numeric scalar
#' @param a numeric scalar
#' @param jmax integer scalar which regulates the accuracy of the result
#' @param cut.point scalar number which regulates the behaviour of the algorithm
#' @return A number between 0 and 1.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @examples
# # OwenT(h,a) = OwenT(-h,a)
#' OwenT(2,1) == OwenT(-2,1)
# # OwenT(0,a) = atan(a)/2pi
#' a <- runif(1, -1000, 1000)
#' OwenT(0,a) - atan(a)/(2*pi)
# # OwenT(h,1) = Phi(h)(1-Phi(h))/2
#' h <- runif(1, -3, 3)
#' OwenT(h,1) - pnorm(h)*(1-pnorm(h))/2
# # OwenT(h,Inf) = (1-Phi(|h|))/2 :
#' OwenT(1,10000) - (1-pnorm(abs(1)))/2
#' OwenT(1,Inf) == (1-pnorm(abs(1)))/2
#' @export
OwenT <- function (h, a, jmax = 50L, cut.point = 8)
{
  if(isNotPositiveInteger(jmax)){
    stop("`jmax` must be an integer >=1.")
  }
  if(cut.point <= 0){
    stop("`cut.point` must be a strictly positive number")
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
  return(RcppOwenT(h, a, jmax, cut.point))
}
