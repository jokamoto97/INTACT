#' Multi-INTACT EM algorithm log-likelihood function.
#'
#' @param w A vector of prior parameter estimates.
#' @param bf A vector of Bayes factors.
#' @param fp_coloc A vector of aggregated colocalization evidence. The default
#' is the max of all GLCPs, truncated at 0.05.
#' @return A scalar log likelihood.



.bf_loglik_coloc<-function(w,bf,fp_coloc){
  K<-length(w)
  n<-length(bf)/K
  loglik<-0
  sumloglik<-0
  for (i in 1:n){
    bfs <- bf[((i-1)*K+1):(i*K)]
    med <- c(exp(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i]) + bfs[1] -
                   .bf_weighted_sum_coloc(w,bf,fp_coloc,i)),
             exp(log(w[2:K]*fp_coloc[i]) + bfs[2:K] -
                   .bf_weighted_sum_coloc(w,bf,fp_coloc,i)))
    loglik <- sum(c(log(w[1]*fp_coloc[i]+ 1 - fp_coloc[i])*med[1],
                    log(w[2:K])*med[2:K]))
    sumloglik<-sumloglik+loglik
  }
  return(sumloglik)
}
