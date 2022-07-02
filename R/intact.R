#' Compute the posterior probability that a gene may be causal, given a gene's TWAS scan z-score and colocalization probability.
#'
#' @param GLCP_vec A vector of colocalization probabilities for the genes of interest
#' @param prior_fun A function to transform a colocalization probability into a prior.
#' Options are linear, step, expit, and hybrid.
#' @param z_vec A vector of TWAS scan z-scores for the genes of interest. The order of
#' genes must match GLCP_vec.
#' @param t A hard threshold for the GLCP. Values below this number will be shrunk to zero.
#'  This argument is used in the user-specified prior function.
#' @param K A vector of values over which Bayesian model averaging is performed.
#' @return The vector of posteriors.
#' @export
#' @examples
#' intact(GLCP_vec=simdat$GLCP, prior_fun=linear, z_vec = simdat$TWAS_z,t = 0.05)





intact <- function(GLCP_vec, prior_fun, z_vec, t=NULL, K = c(1,2,4,8,16)){

  #Compute log10(BF) grid from TWAS z score for each K value

  bf_grid_log10 <- sapply(K,FUN = function(K,z_vec){
    out <- 0.5*log10(1/(1+K)) + (0.5*(z_vec^2) * (K/(1+K)))*log10(exp(1))
    return(out)
  },
  z_vec=z_vec)


  #Perform Bayesian model averaging using log sum exp trick

  bf_log10 <- apply(bf_grid_log10,FUN=function(bf_vec){
    x_star <- max(bf_vec)
    logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
    out <- logsumexp - log10(length(bf_vec))
    return(out)
  },
  MARGIN = 1)

  #Compute prior from GLCPs

  prior <- prior_fun(GLCP_vec,t,u = pi1_fun(z_vec = z_vec,lambda = 0.5))

  #compute log10(posterior odds)

  log10odds <- bf_log10 + log10((prior/(1-prior)))

  #compute posterior probability

  post <- 1/ (1 + 10^(-log10odds))

  return(post)
}
