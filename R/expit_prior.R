#' Transform a gene colocalization probability (GLCP) to a prior to be used in the evidence
#' integration procedure. There are four prior function options, including expit, linear,
#' step, and expit-linear hybrid.
#'
#' @param GLCP A gene colocalization probability
#' @param t A hard threshold for the GLCP. Values below this number will be shrunk to zero.
#' Default is 0.05.
#' @param D A curvature shrinkage parameter. Lower values of D will result in a steeper curve.
#' Default is 0.1
#' @param u A factor between 0 and 1 by which the prior function is scaled.
#' @return The value of the prior.
#' @export
#' @examples
#' expit(0.2, 0.05, 1)



expit <- function(GLCP,t=0.05,D = 0.1,u=1){
  out <- ifelse(GLCP < t, 0,
                1/(1 + exp(-1*(GLCP-0.5)/D)))
  out <- out * u
  return(out)
}

