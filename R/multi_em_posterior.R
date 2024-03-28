#' Compute gene product relevance probabilities using prior parameter estimates
#' and Bayes factors.
#'
#' @param w A vector of prior parameter estimates.
#' @param bf A vector of Bayes factors.
#' @param fp_coloc A vector of aggregated colocalization evidence. The default
#' is the max of all GLCPs, truncated at 0.05.
#' @return A vector of posterior probabilities.



.multi_em_posteriors<-function(w,bf,fp_coloc){
  K<-length(w)
  wnew<-rep(NA,K)
  n<-length(bf)/K
  ei<-matrix(NA,n,K)
  for (i in 1:n){
    bfs <- bf[((i-1)*K+1):(i*K)]
    ei[i,] <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] -
                      .bf_weighted_sum_coloc(w,bf,fp_coloc,i)),
                exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] -
                      .bf_weighted_sum_coloc(w,bf,fp_coloc,i)))
  }
  return(ei)
}

