#' Transform a gene colocalization probability (GLCP) to a prior to be used in the
#' evidence integration procedure. There are four prior function options, including expit,
#' linear, step, and expit-linear hybrid.
#'
#' @param GLCP A gene colocalization probability
#' @param t A hard threshold for the GLCP. Values below this number will be shrunk to zero.
#' Default is 0.05.
#' @param u A factor between 0 and 1 by which the prior function is scaled.
#' @return The value of the prior.
#' @export
#' @examples
#' linear(0.2, 0.05, 1)
#' linear(c(0.01,0.2,0.9))




linear <- function(GLCP, t=0.05,u=1){
  out <- ifelse(GLCP > t, GLCP, 0)
  out <- out * u
  return(out)
}

