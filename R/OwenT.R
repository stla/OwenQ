#' @title Owen T-function
#' @description Evaluates the Owen T-function.
#' @param h numeric scalar
#' @param a numeric scalar
#' @return A number between 0 and 1.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib OwenQ
#' @export
#' @details If \eqn{0\lea\le1}, and \eqn{0\leh\lec}, where \eqn{c} is
#' the cut point, a series expansion is used.
#' It is truncated after the \code{jmax}-th term.
#' If \eqn{0\lea\le1}, and \eqn{h>c}, an asymptotic approximation is used.
#' Otherwise, the properties of the Owen T-function are exploited
#' to come down to the case \eqn{0\lea\le1}.
#' See the reference for more information.
#' @references
#' Owen, D. B. (1956).
#' Tables for computing bivariate normal probabilities.
#' \emph{Ann. Math. Statist.} \bold{27}, 1075-1090.
#' @examples
#' integrate(function(x) pnorm(1+2*x)^2*dnorm(x), lower=-Inf, upper=Inf)
#' pnorm(1/sqrt(5)) - 2*OwenT(1/sqrt(5), 1/3)
OwenT <- function (h, a)
{
  if (!is.numeric(a) || length(a) > 1L)
    stop("`a` must be a scalar number")
  if (!is.numeric(h) || length(h) > 1L)
    stop("`h` must be a scalar number")
  if (is.na(a))
    stop("parameter `a` is NA")
  return(RcppOwenT(h, a))
}
