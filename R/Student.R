Ssequences <- function(n, a, b, d){
  A <- M <- numeric(n)
  sb <- sqrt(b)
  M[1L] <- a * sb * dnorm(d*sb) * pnorm(d*a*sb)
  A[2L] <- 1
  M[2L] <- b * (d * a * M[1L] + a * dnorm(d) / sqrt(2*pi))
  for(k in 3L:n){
    A[k] <- 1 /(k-2)/A[k-1L]
    M[k] <- (k-2)/(k-1) * b * (A[k-1L] * d * a * M[k-1L] + M[k-2L])
  }
  return(M)
}

#' @title Student CDF with integer number of degrees of freedom
#' @description Cumulative distribution function of the noncentrel Student
#' distribution with an integer number of degrees of freedom.
#' @param q quantile
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param delta noncentrality parameter
#' @param jmax,cut.point passed to \code{\link{OwenT}} (when \code{nu} is odd)
#' @return Numeric value, the CDF evaluated at \code{q}.
#' @export
#' @importFrom stats pnorm dnorm
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' ptOwen(2, 3) - pt(2, 3)
#' ptOwen(2, 3, delta=1) - pt(2, 3, ncp=1)
ptOwen <- function(q, nu, delta=0, jmax=50L, cut.point=8){
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(q) || is.infinite(delta)){
    stop("Parameters must be finite.")
  }
  a <- sign(q)*sqrt(q*q/nu)
  b <- nu/(nu+q*q)
  nu <- as.integer(nu)
  if(nu==2L){
    sB <- sqrt(b)
    asB <- sign(q)*sqrt(q*q/(nu+q*q))
    return(pnorm(-delta) + sqrt(2*pi) *
             (asB * dnorm(delta*sB) * pnorm(delta*asB)))
  }
  if(nu%%2L==1L){
    sB <- sqrt(b)
    C <- pnorm(-delta*sB) + 2*.OwenT(delta*sB,a, jmax=jmax, cut.point=cut.point)
    if(nu==1L){
      return(C)
    }
    if(nu==3L){
      ab <- a*b
      asB <- sign(q)*sqrt(q*q/(nu+q*q))
      return(C + 2 * ab * (delta * asB * dnorm(delta*sB) * pnorm(delta*asB)
                     + exp(-delta*delta/2)/(2*pi)))
    }
    return(C + 2*sum(Ssequences(nu-1, a, b, delta)[seq(2L, nu-1L, by=2L)]))
  }
  return(pnorm(-delta) + sqrt(2*pi) * sum(Ssequences(nu-1, a, b, delta)[seq(1L, nu-1L, by=2L)]))
}

