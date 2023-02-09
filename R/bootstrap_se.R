#' Compute bootstrap standard errors for alpha MLEs.
#'
#' @param d_vec A vector of gene set annotations for the genes of interest.
#' Entries should be integer(1) if the gene is annotated and integer(0)
#' otherwise.
#' @param pprobs A vector of posterior probabilities for each gene estimated
#' from
#' the intact function. Gene order should match d_vec.
#' @param reps Number of bootstrap samples.
#' @return MLEs for alpha0 and alpha1 from bootstrap samples.



.enrich_bootstrap_se <- function(pprobs, d_vec, reps = 100){

  #generate bootstrap samples indices

  samples <- replicate(n = 100,
                       expr = sample(seq(length(pprobs)),
                                     size = length(pprobs), replace=TRUE))

  #compute gene set enrichment for each bootstrap sample

  res <- apply(samples, 2, FUN=function(x){
    boot_dat <- cbind(pprobs,d_vec)[x,]
    return(.em_est(pprobs=boot_dat[,1], d_vec=boot_dat[,2]))
  })
}
