#' Compute Multi-INTACT prior parameter estimates and gene product relevance
#' probabilities.
#'
#' @param df A data frame with marginal z-scores for the TWAS analyses of
#' gene products 1 and 2, as well as the multivariate Wald chi-square statistic
#' from the joint analysis.
#' @param pi_init Initialization of prior parameters to be estimated. The first
#' entry should be the qvalue estimate for pi0. The second parameter should
#' be the gene-product-1-only parameter; the third should be the
#' gene-product-2-only parameter, and the fourth should be the
#' gene-product-1+2-parameter.
#' @param chisq_vec A numeric vector of multivariate Wald chi-square test
#' statistics. The order of this vector should match z_1, z_2, and fp_coloc.
#' @param chisq_dof A numeric vector containing the chi-square test statistic
#' degrees of freedom under the null. Default is 2.
#' @param z_1 A numeric vector of TWAS test statistics for gene product 1.
#' The order of this vector should match chisq_vec, z_2, and fp_coloc.
#' @param z_2 A numeric vector of TWAS test statistics for gene product 2.
#' The order of this vector should match chisq_vec, z_1, and fp_coloc.
#' @param fp_coloc A vector of aggregated colocalization evidence. The default
#' is the max of all GLCPs, truncated at 0.05. The order should match
#' chisq_vec, z_1, and z_2.
#' @return A list containing a data frame with the model posteriors
#' (posterior_1, posterior_2, and posterior_12),gene product relevance
#' probabilities (GPRP_1 and GPRP_2), a vector of prior parameter estimates and
#' a Boolean indicating convergence of the EM algorithm.
#' @importFrom tidyr pivot_longer



.multi_prior_estimation <- function(df,
                                    pi_init,
                                    chisq_vec,
                                    chisq_dof = 2,
                                    z_1,
                                    z_2,
                                    fp_coloc){

  df$xwas_z <- qnorm(pchisq(chisq_vec,df = chisq_dof,
                            lower.tail = FALSE)/2,lower.tail = FALSE)

  # Compute Bayes factors for each model

  df$BF_12 <- .wakefield_bf_z_ln(df$xwas_z)
  df$BF_1 <- .wakefield_bf_z_ln(df$z_1)
  df$BF_2 <- .wakefield_bf_z_ln(df$z_2)
  df$BF_0 <- log(1)
  bf_df <- pivot_longer(df,cols = c(BF_0,BF_1,BF_2,BF_12),
                        names_to = "BF_Type",
                        values_to = "BF")

  #Estimate prior parameters

  pi_hat_coloc_pi0 <- SQUAREM::squarem(par = pi_init,
                                       fixptfn = .bf_em_coloc_pi0,
                                       objfn = .bf_loglik_coloc,
                                       control = list(tol = 1.e-08,
                                                      minimize=FALSE,
                                                      maxiter=50),
                                       bf = bf_df$BF,
                                       fp_coloc = fp_coloc)

  pi_hats <- round(pi_hat_coloc_pi0$par,digits = 4)

  converged <- pi_hat_coloc_pi0$convergence

  #If EM algorithm failed to converge, use pi_init to compute posteriors

  if(converged == FALSE){

    pi_hats <- pi_init

  }

  #Compute posteriors for each model

  df$posterior_0 <- .multi_em_posteriors(w = pi_hats,
                                         bf = bf_df$BF,
                                         fp_coloc = fp_coloc)[,1]
  df$posterior_1 <- .multi_em_posteriors(w = pi_hats,
                                         bf = bf_df$BF,
                                         fp_coloc = fp_coloc)[,2]
  df$posterior_2 <- .multi_em_posteriors(w = pi_hats,
                                         bf = bf_df$BF,
                                         fp_coloc = fp_coloc)[,3]
  df$posterior_12 <- .multi_em_posteriors(w = pi_hats,
                                         bf = bf_df$BF,
                                         fp_coloc = fp_coloc)[,4]

  df$GPRP_1 <- df$posterior_1 + df$posterior_12

  df$GPRP_2 <- df$posterior_2 + df$posterior_12

  #Return list of results

  out <- list(df,pi_hats,converged)

  return(out)
}

