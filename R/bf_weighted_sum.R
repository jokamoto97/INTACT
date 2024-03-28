#' Helper function for EM algorithm to compute weighted sum of Bayes factors.
#'
#' @param w A vector of prior parameter estimates.
#' @param bf A vector of Bayes factors.
#' @param fp_coloc A vector of aggregated colocalization evidence. The default
#' is the max of all GLCPs, truncated at 0.05.
#' @param i The index for the ith gene.
#' @return A scalar weighted sum of Bayes factors



.bf_weighted_sum_coloc<-function(w,bf,fp_coloc,i){
  K <- length(w)
  fp_coloc <- fp_coloc[i]
  bf.sum <- 0
  bf.gene <- bf[((i-1)*K+1):(i*K)]
  bf.m <- max(bf.gene)
  bf.null <- bf.gene[1]
  bf.alt <- bf.gene[2:length(bf.gene)]
  bf.sum <- sum(c((w[1]*fp_coloc + 1 - fp_coloc)*exp(bf.null-bf.m),
                  w[2:K]*fp_coloc*exp(bf.alt-bf.m)))
  return(bf.m+log(bf.sum))
}
