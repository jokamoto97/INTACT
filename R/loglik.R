#' A log likelihood function for the expectation-maximization algorithm.
#' Used as an argument for objfn in the squarem function.
#'
#' @param d_vec A vector of gene set annotations for the genes of interest.
#' Entries should be integer(1) if the gene is annotated and integer(0)
#' otherwise.
#' @param pprobs A vector of posterior probabilities for each gene estimated
#' from the
#' intact function. Gene order should match d_vec.
#' @param alpha A vector containing the current estimates of the enrichment
#' parameters
#' alpha0 and alpha1.
#' @return Log likelihood evaluated at the current estimates of alpha0 and
#' alpha1.
#' @examples
#' data(simdat)
#' .logistic_loglik(d_vec = sample(c(0,1),1197,replace=TRUE),
#' pprobs = intact(GLCP_vec=simdat$GLCP, prior_fun=linear,
#' z_vec = simdat$TWAS_z, t = 0.05),
#' alpha = c(0,0))



.logistic_loglik <- function(alpha,d_vec,pprobs){

  #Empirical Bayes estimate for prior
  pi <- mean(pprobs)

  #account for genes with pprob == 1 so that likelihood does not become Inf
  pprobs <- ifelse(pprobs > 0, pprobs - 1.e-8, pprobs)

  BF <- (pprobs/(1-pprobs))/(pi/(1-pi))

  liklhd <- (1/(1+exp(-1*(alpha[1] + alpha[2]*d_vec))))*BF +
    1/(1+exp((alpha[1] + alpha[2]*d_vec)))

  loglik <- log(liklhd)

  return(sum(loglik))
}

