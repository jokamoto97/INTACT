#' Estimate pi1 from TWAS scan z-scores.
#'
#' @param z_vec A vector of TWAS scan z-scores.
#' @param lambda A value between 0 and 1. The density of TWAS scan z-scores
#' should be
#' flat at lambda. Set to 0.5 as default.
#' @return A scalar estimate for pi1.
#' @export
#' @examples
#' data(simdat)
#' pi1_fun(simdat$TWAS_z)
#' @importFrom stats pnorm qnorm



pi1_fun <- function(z_vec,lambda = 0.5){

  p_vec <- 2*pnorm(abs(z_vec),lower.tail = FALSE)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)
}
