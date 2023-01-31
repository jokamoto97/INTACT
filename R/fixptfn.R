#' A fixed-point mapping for the expectation-maximization algorithm.
#' Used as an argument for fixptfn in the squarem function.
#'
#' @param d_vec A vector of gene set annotations for the genes of interest.
#' Entries should be integer(1) if the gene is annotated and integer(0)
#' otherwise.
#' @param pprobs A vector of posterior probabilities for each gene estimated
#' from the
#' intact function. Gene order should match d_vec.
#' @param alpha A vector containing the current estimates of the enrichment
#' parameters
#' alpha0 and alpha1$.
#' @return Updated estimates of alpha0 and alpha1.
#' @export
#' @examples
#' data(simdat)
#' logistic_em(d_vec = sample(c(0,1),1197,replace=TRUE),
#' pprobs = intact(GLCP_vec=simdat$GLCP, prior_fun=linear,
#' z_vec = simdat$TWAS_z, t = 0.05),
#' alpha = c(0,0))

logistic_em <- function(d_vec,pprobs,alpha){

  #new alpha vector

  alpha_new <- rep(NA,2)

  pi <- mean(pprobs)

  #Induced Bayes factor vector (posterior odds/prior odds), on log10 scale to prevent overflow:

  log10BF <- log10(pprobs/(1-pprobs)) -log10(pi/(1-pi))

  ####E-step####
  pi_t <- 1/(1 + exp(-1*(alpha[1] + alpha[2]*d_vec)))

  #compute log10(posterior odds)

  log10odds <- log10BF + log10((pi_t/(1-pi_t)))

  #Update gamma estimate

  Egamma <- 1 / (1 + 10^(-log10odds))

  #M step

  C00 <- sum((1-Egamma)*(1-d_vec)) + 1/length(d_vec)

  C10 <- sum(Egamma*(1-d_vec)) + 1/length(d_vec)

  C01 <- sum((1-Egamma)*d_vec) + 1/length(d_vec)

  C11 <- sum(Egamma * d_vec) + 1/length(d_vec)

  #Update alpha

  alpha0 <- log(C10/C00)

  alpha1 <- log(C00*C11/(C10*C01))

  alpha_new <- c(alpha0,alpha1)

  diff <- sqrt(crossprod(alpha_new-alpha))

  alpha <- alpha_new

  return(alpha_new)
}
