#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t number, positive or negative, possibly infinite
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (upper bound of the integral) vector of finite positive numbers, with the same length as \code{delta}
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @importFrom stats pgamma
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(nu=5, t=3, delta=2, R=100)
#' pt(q=3, df=5, ncp=2)
OwenQ1 <- function(nu, t, delta, R){
  J <- length(delta)
  if(J != length(R)){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0 || R==Inf)){
    stop("`R` must be a finite positive number.")
  }
  if(nu==Inf){
    return(numeric(J))
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(delta))){
    out <- numeric(J)
    if(any(minf <- delta==-Inf)){
      wminf <- which(minf)
      out[wminf] <- pgamma(R[wminf]^2/2, nu/2, lower.tail=TRUE)
    }
    if(!all(inf <- is.infinite(delta))){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenQ1(nu, t, delta[noninf], R[noninf])
    }
    return(out)
  }
  RcppOwenQ1(nu, t, delta, R)
}

#' @title Second Owen Q-function
#' @description Evaluates the second Owen Q-function (integral from \eqn{R} to \eqn{\infty})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t number, positive or negative, possibly infinite
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (lower bound of the integral) vector of finite positive numbers,
#' with the same length as \code{delta}
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral
#' from \eqn{R} to \eqn{\infty}.
#' @export
#' @importFrom stats pgamma pnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' # OwenQ1(nu, t, delta, R) + OwenQ2(nu, t, delta, R) = pt(t, nu, delta)
#' OwenQ1(nu=5, t=3, delta=2, R=1) + OwenQ2(nu=5, t=3, delta=2, R=1)
#' pt(q=3, df=5, ncp=2)
OwenQ2 <- function(nu, t, delta, R){
  J <- length(R)
  if(length(delta) != J){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0 || R==Inf)){
    stop("`R` must be a finite positive number.")
  }
  if(nu==Inf){
    return(pnorm(t, mean=delta))
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(delta))){
    out <- numeric(J)
    if(any(minf <- delta==-Inf)){
      wminf <- which(minf)
      out[wminf] <- pgamma(R[wminf]^2/2, nu/2, lower.tail=FALSE)
    }
    if(!all(inf <- is.infinite(delta))){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenQ2(nu, t, delta[noninf], R[noninf])
    }
    return(out)
  }
  RcppOwenQ2(nu, t, delta, R)
}
