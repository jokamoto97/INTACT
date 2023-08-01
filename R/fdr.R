#' Bayesian FDR control for INTACT output
#'
#' @param posterior A vector of posterior probabilities for each gene estimated
#' from the
#' intact function.
#' @param alpha A numeric target FDR control level.
#' @return An n x 2 data frame where the first column is the inputted
#' posterior probabilities, and the second is a Boolean vector denoting
#' significance at the specified target control level.
#' @export
#' @examples
#' data(simdat)
#' fdr_rst(simdat$GLCP)


fdr_rst<-function(posterior, alpha=0.05){

  gene_num <- seq(1,length(posterior))

  lfdr_sort = sort(1-posterior)

  FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))

  thresh = 1 - lfdr_sort[max(which(FDR<=alpha))]

  rej_gene = as.numeric(gene_num[which(posterior>=thresh)])

  out_tmp <- rep(FALSE,length(posterior))

  out_tmp[rej_gene] <- TRUE

  out <- data.frame("posterior" = posterior,
                    "sig" = out_tmp)

  return(out)
}
