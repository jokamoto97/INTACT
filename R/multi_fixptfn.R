#' Multi-INTACT EM algorithm fixed-point function.
#'
#' @param w A vector of prior parameter estimates.
#' @param bf A vector of Bayes factors.
#' @param fp_coloc A vector of aggregated colocalization evidence. The default
#' is the max of all GLCPs, truncated at 0.05.
#' @return A vector of updated parameter estimates.



.bf_em_coloc_pi0<-function(w,bf,fp_coloc){
  K<-length(w)
  wnew<-rep(NA,K)
  n<-length(bf)/K
  #E-step
  ei<-matrix(NA,n,K)
  for (i in 1:n){
    bfs <- bf[((i-1)*K+1):(i*K)]
    ei[i,] <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] -
                      .bf_weighted_sum_coloc(w,bf,fp_coloc,i)),
                exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] -
                      .bf_weighted_sum_coloc(w,bf,fp_coloc,i)))
  }
  #M-step
  wnew[1] <- w[1]
  diff <- 1 - wnew[1]
  tmp <- colSums(ei[,2:K],na.rm = TRUE)
  scaling <- sum(tmp)/diff
  wnew[2:K] <- tmp/scaling
  return(wnew)
}
