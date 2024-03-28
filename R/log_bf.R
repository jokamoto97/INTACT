#' A function to compute log Bayes factors from z-statistics using the
#' Wakefield formula
#'
#' @param z_vec A vector of TWAS scan z-scores for the genes of interest.
#' @param K A vector of values over which Bayesian model averaging is performed.
#' @return log Bayes factors
#'


.wakefield_bf_z_ln <- function(z_vec, K = c(1,2,4,8,16)){

  bf_grid_log <- vapply(X=matrix(K),FUN.VALUE=numeric(length(z_vec)),
                        FUN = function(K,z_vec){
                          out <- 0.5*log(1/(1+K)) + (0.5*(z_vec^2) * (K/(1+K)))
                          return(out)
                        },
                        z_vec=z_vec)


  #Perform Bayesian model averaging using log sum exp trick

  bf_log <- apply(bf_grid_log,FUN=function(bf_vec){
    x_star <- max(bf_vec)
    logsumexp <- x_star + log(sum(exp(bf_vec-x_star)))
    out <- logsumexp - log(length(bf_vec))
    return(out)
  },
  MARGIN = 1)

  return(bf_log)
}

