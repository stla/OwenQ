#' @title Owen distribution functions when \eqn{\delta_1>\delta_2}
#' @description Evaluates the Owen cumulative distribution function
#' when the noncentrality parameters satisfy \eqn{\delta_1>\delta_2} and
#' the number of degrees of freedom is integer.
#' \itemize{
#' \item \code{powen1} evaluates \eqn{P(T_1\le t_1, T_2 \le t_2)}
#' (Owen's equality 8)
#' \item \code{powen2} evaluates \eqn{P(T_1\le t_1, T_2 \ge t_2)}
#' (Owen's equality 9)
#' \item \code{powen3} evaluates \eqn{P(T_1\ge t_1, T_2 \ge t_2)}
#' (Owen's equality 10)
#' \item \code{powen4} evaluates \eqn{P(T_1\ge t_1, T_2 \le t_2)}
#' (Owen's equality 11)
#' }
#' @name powen
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom;
#' infinite allowed
#' @param t1,t2 two finite numbers, positive or negative
#' @param delta1,delta2 two vectors of possibly infinite numbers with the same length,
#' the noncentrality parameters;
#' must satisfy \code{delta1>delta2}
#' @param jmax,cutpoint parameters controlling the algorithm for the Owen-T function;
#' see \code{\link{OwenT}} (used only when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}.
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @seealso Use \code{\link{OwenCDF}} for general values of \code{delta1} and \code{delta2}.
#' @references
#' Owen, D. B. (1965).
#' A special case of a bivariate noncentral t-distribution.
#' \emph{Biometrika} \bold{52}, 437-446.
#' @examples
#' nu=5; t1=2; t2=1; delta1=3; delta2=2
#' # Wolfram integration gives 0.1394458271284726
#' ( p1 <- powen1(nu, t1, t2, delta1, delta2) )
#' # Wolfram integration gives 0.0353568969628651
#' ( p2 <- powen2(nu, t1, t2, delta1, delta2) )
#' # Wolfram integration gives 0.806507459306199
#' ( p3 <- powen3(nu, t1, t2, delta1, delta2) )
#' # Wolfram integration gives 0.018689824158
#' ( p4 <- powen4(nu, t1, t2, delta1, delta2) )
#' # the sum should be 1
#' p1+p2+p3+p4
NULL

#' @rdname powen
#' @export
powen1 <- function(nu, t1, t2, delta1, delta2, jmax=50L, cutpoint=8){
  J <- length(delta1)
  if(J != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2 & is.finite(delta1) & is.finite(delta2))){
    stop("`delta1` must be >`delta2`.")
  }
  if(is.infinite(t1) || is.infinite(t2)){
    stop("`t1` and `t2` must be finite.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(t1<=t2){
    return(ptOwen(t1, nu, delta1))
  }
  if(nu == Inf){
    return(pnorm(t1, mean=delta1) - pmax(0, pnorm(t1, mean=delta1)-pnorm(t2, mean=delta2)))
  }
  if(any(inf <- (is.infinite(delta1) | is.infinite(delta2)))){
    out <- numeric(J)
    if(any(inf2 <- delta2==-Inf)){
      winf2 <- which(inf2)
      out[winf2] <- ptOwen(t1, nu, delta1[winf2])
    }
    if(!all(inf)){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenCDF1(nu, t1, t2, delta1[noninf], delta2[noninf], jmax, cutpoint)
    }
    return(out)
  }
  RcppOwenCDF1(nu, t1, t2, delta1, delta2, jmax, cutpoint)
}

#' @rdname powen
#' @export
powen2 <- function(nu, t1, t2, delta1, delta2, jmax=50L, cutpoint=8){
  J <- length(delta1)
  if(J != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2 & is.finite(delta1) & is.finite(delta2))){
    stop("`delta1` must be >`delta2`.")
  }
  if(is.infinite(t1) || is.infinite(t2)){
    stop("`t1` and `t2` must be finite.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(t1<=t2){
    return(0)
  }
  if(nu == Inf){
    #return(-pnorm(t2, mean=delta2)+pnorm(t1, mean=delta1) + pmax(0, pnorm(t2, mean=delta2)-pnorm(t1, mean=delta1)))
    return(pmax(0, pnorm(t1, mean=delta1)-pnorm(t2, mean=delta2)))
  }
  if(any(inf <- (is.infinite(delta1) | is.infinite(delta2)))){
    out <- numeric(J)
    if(!all(inf)){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenCDF2(nu, t1, t2, delta1[noninf], delta2[noninf], jmax, cutpoint)
    }
    return(out)
  }
  RcppOwenCDF2(nu, t1, t2, delta1, delta2, jmax, cutpoint)
}

#' @rdname powen
#' @export
powen3 <- function(nu, t1, t2, delta1, delta2, jmax=50L, cutpoint=8){
  J <- length(delta1)
  if(J != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2 & is.finite(delta1) & is.finite(delta2))){
    stop("`delta1` must be >`delta2`.")
  }
  if(is.infinite(t1) || is.infinite(t2)){
    stop("`t1` and `t2` must be finite.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  # if(any(lower <- (delta1<delta2))){
  #   out <- numeric(J)
  #   out[lower] <- 1 - ptOwen(t1, nu, delta1, jmax, cutpoint)
  #   if(length(nonlower <- which(!lower))){
  #     out[nonlower] <- RcppOwenCDF3(nu, t1, t2, delta1[nonlower], delta2[nonlower],
  #                                   jmax, cutpoint)
  #   }
  #   return(out)
  # }
  if(t1<=t2){
    return(1-ptOwen(t2, nu, delta2, jmax, cutpoint))
  }
  if(nu == Inf){ # to simplify ?
    return(1-pnorm(t1, mean=delta1) - pmax(0, pnorm(t2, mean=delta2)-pnorm(t1, mean=delta1)))
  }
  if(any(inf <- (is.infinite(delta1) | is.infinite(delta2)))){
    out <- numeric(J)
    inf1 <- which(is.infinite(delta1))
    if(length(inf1)){
      out[inf1] <- ifelse(delta1[inf1]==Inf, 1-ptOwen(t2, nu, delta2[inf1], jmax, cutpoint), 0)
    }
    if(!all(inf)){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenCDF3(nu, t1, t2, delta1[noninf], delta2[noninf], jmax, cutpoint)
    }
    return(out)
  }
  RcppOwenCDF3(nu, t1, t2, delta1, delta2, jmax, cutpoint)
}

#' @rdname powen
#' @export
powen4 <- function(nu, t1, t2, delta1, delta2, jmax=50L, cutpoint=8){
  J <- length(delta1)
  if(J != length(delta2)){
    stop("`delta1` and `delta2` must have the same length.")
  }
  if(any(delta1<=delta2 & is.finite(delta1) & is.finite(delta2))){
    stop("`delta1` must be >`delta2`.")
  }
  if(is.infinite(t1) || is.infinite(t2)){
    stop("`t1` and `t2` must be finite.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  # if(any(lower <- (delta1<delta2))){
  #   out <- numeric(J)
  #   #out[lower] <- 0
  #   if(length(nonlower <- which(!lower))){
  #     out[nonlower] <- RcppOwenCDF4(nu, t1, t2, delta1[nonlower], delta2[nonlower],
  #                                   jmax, cutpoint)
  #   }
  #   return(out)
  # }
  if(t1<=t2){
    return(ptOwen(t2, nu, delta2, jmax, cutpoint)-ptOwen(t1, nu, delta1, jmax, cutpoint))
  }
  if(nu == Inf){
    return(pmax(0, pnorm(t2, mean=delta2)-pnorm(t1, mean=delta1)))
  }
  if(any(inf <- (is.infinite(delta1) | is.infinite(delta2)))){
    out <- numeric(J)
    inf1 <- which(delta1==Inf)
    if(length(inf1)){
      out[inf1] <- ptOwen(t2, nu, delta2[inf1], jmax, cutpoint)
    }
    minf2 <- setdiff(which(delta2==-Inf), inf1)
    if(length(minf2)){
      out[minf2] <- 1 - ptOwen(t1, nu, delta1, jmax, cutpoint)
    }
    if(!all(inf)){
      noninf <- which(!inf)
      out[noninf] <- RcppOwenCDF4(nu, t1, t2, delta1[noninf], delta2[noninf], jmax, cutpoint)
    }
    return(out)
  }
  RcppOwenCDF4(nu, t1, t2, delta1, delta2, jmax, cutpoint)
}
