#' @title Owen distribution functions
#' @description Evaluates the Owen cumulative distribution function
#' for an integer number of degrees of freedom.
#' \itemize{
#' \item \code{psbt1} evaluates \eqn{P(T_1\le t_1, T_2 \le t_2)}
#' \item \code{psbt2} evaluates \eqn{P(T_1\le t_1, T_2 \ge t_2)}
#' \item \code{psbt3} evaluates \eqn{P(T_1\ge t_1, T_2 \ge t_2)}
#' \item \code{psbt4} evaluates \eqn{P(T_1\ge t_1, T_2 \le t_2)}
#' }
#' @name psbt
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom;
#' infinite allowed
#' @param t1,t2 two finite numbers, positive or negative
#' @param delta1,delta2 two vectors of possibly infinite numbers with the same length,
#' the noncentrality parameters
#' @return A vector of numbers between \eqn{0} and \eqn{1}.
#' @note When the number of degrees of freedom is odd, the procedure resorts to
#' the Owen T-function (\code{\link{OwenT}}).
#' @seealso It is better to use \code{\link{powen}} if \code{delta1>delta2}.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' nu=5; t1=1; t2=2; delta1=2; delta2=3
#' ( p1 <- psbt1(nu, t1, t2, delta1, delta2) )
#' ( p2 <- psbt2(nu, t1, t2, delta1, delta2) )
#' ( p3 <- psbt3(nu, t1, t2, delta1, delta2) )
#' ( p4 <- psbt4(nu, t1, t2, delta1, delta2) )
#' # the sum should be 1
#' p1+p2+p3+p4
NULL

#' @rdname psbt
#' @export
psbt1 <- function(nu, t1, t2, delta1, delta2){
  L <- length(delta1)
  if(L != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  out <- numeric(L)
  higher <- (delta1 > delta2) | is.infinite(delta1) | is.infinite(delta2)
  if(J <- length(whigher <- which(higher))){
    out[whigher] <- powen1(nu, t1, t2, delta1[whigher], delta2[whigher])
  }
  if(J < L){
    equal <- delta1==delta2
    if(K <- length(wequal <- which(equal))){
      out[wequal] <- ptOwen(min(t1,t2), nu, delta1[wequal])
    }
    if(J+K < L){
      wlower <- which(!(higher|equal))
      if(t1 < t2){
        out[wlower] <- RcppOwenCDF1(nu, t2, t1, delta2[wlower], delta1[wlower])
      }else{
        out[wlower] <- ptOwen(t2, nu, delta2[wlower])
      }
    }
  }
  out
}

#' @rdname psbt
#' @export
psbt2 <- function(nu, t1, t2, delta1, delta2){
  L <- length(delta1)
  if(L != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  out <- numeric(L)
  higher <- (delta1 > delta2) | is.infinite(delta1) | is.infinite(delta2)
  if(J <- length(whigher <- which(higher))){
    out[whigher] <- powen2(nu, t1, t2, delta1[whigher], delta2[whigher])
  }
  if(J < L){
    equal <- delta1==delta2
    if(K <- length(wequal <- which(equal))){
      if(t2 < t1){
        out[wequal] <- ptOwen(t1, nu, delta1[wequal]) -
          ptOwen(t2, nu, delta1[wequal])
      }
    }
    if(J+K < L){
      wlower <- which(!(higher|equal))
      if(t1 < t2){
        out[wlower] <- RcppOwenCDF4(nu, t2, t1, delta2[wlower], delta1[wlower])
      }else{
        out[wlower] <- ptOwen(t1, nu, delta1[wlower]) - ptOwen(t2, nu, delta2[wlower])
      }
    }
  }
  out
}

#' @rdname psbt
#' @export
psbt3 <- function(nu, t1, t2, delta1, delta2){
  L <- length(delta1)
  if(L != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  out <- numeric(L)
  higher <- (delta1 > delta2) | is.infinite(delta1) | is.infinite(delta2)
  if(J <- length(whigher <- which(higher))){
    out[whigher] <- powen3(nu, t1, t2, delta1[whigher], delta2[whigher])
  }
  if(J < L){
    equal <- delta1==delta2
    if(K <- length(wequal <- which(equal))){
      out[wequal] <- 1 - ptOwen(max(t1,t2), nu, delta1[wequal])
    }
    if(J+K < L){
      wlower <- which(!(higher|equal))
      if(t1 < t2){
        out[wlower] <- RcppOwenCDF3(nu, t2, t1, delta2[wlower], delta1[wlower])
      }else{
        out[wlower] <- 1-ptOwen(t1, nu, delta1[wlower])
      }
    }
  }
  out
}

#' @rdname psbt
#' @export
psbt4 <- function(nu, t1, t2, delta1, delta2){
  L <- length(delta1)
  if(L != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  out <- numeric(L)
  higher <- (delta1 > delta2) | is.infinite(delta1) | is.infinite(delta2)
  if(J <- length(whigher <- which(higher))){
    out[whigher] <- powen4(nu, t1, t2, delta1[whigher], delta2[whigher])
  }
  if(J < L){
    equal <- delta1==delta2
    if(K <- length(wequal <- which(equal))){
      if(t2 > t1){
        out[wequal] <- ptOwen(t2, nu, delta1[wequal]) -
          ptOwen(t1, nu, delta1[wequal])
      }
    }
    if(J+K < L){
      wlower <- which(!(higher|equal))
      if(t1 < t2){
        out[wlower] <- RcppOwenCDF2(nu, t2, t1, delta2[wlower], delta1[wlower])
      }
    }
  }
  out
}
