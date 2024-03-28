#' Compute a gene-level multivariate Wald chi-square statistic using
#' summary-level genetic association and LD data.
#'
#' @param z_vec A 2-vector for which the first entry is the TWAS z-score for
#' the first gene product, and the second entry is the TWAS z-score for the
#' second gene product (e.g. TWAS and PWAS z-scores).
#' @param w A 2 x p matrix of TWAS prediction model weights that are used to
#' generate the statistics required for the z_vec argument. The order of the
#' columns must match the order of the z_vec statistics. Weights should be
#' generated from standardized genotypes.
#' @param R A p x p LD correlation matrix containing pairwise correlations for
#' each SNP included in w.
#' @return The value of approximated multivariate chi-square statistic.
#' @export
#' @examples
#' data(z_sumstats)
#' data(protwt_sumstats)
#' data(exprwt_sumstats)
#' data(ld_sumstats)
#' chisq_sumstat(z_vec = z_sumstats,w = cbind(protwt_sumstats,exprwt_sumstats),
#' R = ld_sumstats)




chisq_sumstat <- function(z_vec, w, R){

  if (!(is.numeric(z_vec))){

    stop("z_vec must be numeric.")

  }
  if(!(is.numeric(w))){

    stop("w must be numeric.")

  }
  if(!(is.numeric(R))){

    stop("R must be numeric.")

  }
  if(ncol(w) != length(z_vec)){

    stop("ncol(w) must match length(z).")

  }
  if(nrow(R) != ncol(R)){

    stop("R must be a square matrix.")

  }
  if(nrow(w) != nrow(R)){

    stop("nrow(w) must match nrow(R).")

  }
  offdiag <- (t(w[,1]) %*% R %*% w[,2])/sqrt(t(w[,1]) %*%
                                              R %*% w[,1] * t(w[,2]) %*%
                                              R %*% w[,2])
  #Predicted gene product correlation matrix

  cor_mat <- cbind(c(1,offdiag),c(offdiag,1))

  inv_cor_mat <- try(solve(cor_mat))

  if("try-error" %in% class(inv_cor_mat)){

    stop("Correlation matrix is not invertible. If the predicted gene products
         are higly correlated, we recommend using z^2 with 1 degree of freedom
         as the chi-square statistic.")

  }

  X2 <- t(z_vec) %*% inv_cor_mat %*% z_vec

  return(X2)
}


