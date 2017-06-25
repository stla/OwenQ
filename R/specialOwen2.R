#' @title Special case of second Owen distribution function
#' @description Evaluation of the second Owen distribution function in a
#' special case (see details).
#' @param nu positive integer, possibly infinite
#' @param t positive number
#' @param delta vector of positive numbers
#'
#' @return A vector of numbers between 0 and 1.
#' @importFrom stats pnorm
#' @export
#' @details The value of \code{spowen2(nu, t, delta)} is the same as the value
#' of \code{powen2(nu, t, -t, delta, -delta)}, but it is evaluated more efficiently.
#' @seealso \code{\link{powen2}}
#' @examples
#' spowen2(4, 1, 2) == powen2(4, 1, -1, 2, -2)
spowen2 <- function(nu, t, delta){
  J <- length(delta)
  if(is.infinite(t)){
    stop("`t` must be finite.")
  }
  if(t<0){
    stop("`t` must be positive")
  }
  if(any(delta<0)){
    stop("`delta` must be positive")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(nu == Inf){
    return(pmax(0, 2*pnorm(t, mean=delta)-1))
  }
  if(any(inf <- (delta==Inf))){
    out <- numeric(J)
    if(!all(inf)){
      noninf <- which(!inf)
      out[noninf] <- RcppSpecialOwenCDF2(nu, t, delta[noninf])
    }
    return(out)
  }
  RcppSpecialOwenCDF2(nu, t, delta)
}
