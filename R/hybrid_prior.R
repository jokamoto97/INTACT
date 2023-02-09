#' Transform a gene colocalization probability (GLCP) to a prior to be used in
#' the
#' evidence integration procedure. There are four prior function options,
#' including expit,
#' linear, step, and expit-linear hybrid.
#'
#' @param GLCP A gene colocalization probability
#' @param t A hard threshold for the GLCP. Values below this number will be
#' shrunk to zero.
#' Default is 0.05.
#' @param D A curvature shrinkage parameter. Lower values of D will result in a
#' steeper curve.
#' Default is 0.1
#' @param u A factor between 0 and 1 by which the prior function is scaled.
#' @param thresholding An option to use hard thresholding or soft thresholding
#' for the prior function. Default is "hard". For soft thresholding, set to
#' "soft".
#' @return The value of the prior.
#' @export
#' @examples
#' hybrid(0.2, 0.05, 1)




hybrid <- function(GLCP, t=0.05, D = 0.1,u=1,thresholding = "hard"){
  if (thresholding == "hard"){
    out <- ifelse(GLCP < t, 0,
                  ifelse(GLCP >= t & GLCP < 0.5,  1/(1 + exp(-1*(GLCP-0.5)/D)),
                        GLCP))
  }
  if (thresholding == "soft"){
    out <- ifelse(GLCP < t, 0,
                  ifelse(GLCP < 0.5 & (1/(1 + exp(-1*(GLCP-0.5)/D)) - t) < 0,0,
                         ifelse(GLCP < 0.5 & (1/(1 + exp(-1*(GLCP-0.5)/D)) - t) > 0,1/(1 + exp(-1*(GLCP-0.5)/D)) - t,
                                ifelse(GLCP >= 0.5 & (GLCP - t) < 0,0,GLCP - t))))
  }
  if (thresholding != "hard" & thresholding != "soft"){
    stop("thresholding must be set to 'hard' or 'soft'.")
  }
  out <- out * u
  return(out)
}


