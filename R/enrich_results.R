#' Compute gene set enrichment estimates with standard errors.
#'
#' @param sig_lev A significance threshold for gene set enrichment hypothesis
#' testing.
#' @param d_vec A vector of gene set annotations for the genes of interest.
#' Entries should be integer(1) if the gene is annotated and integer(0)
#' otherwise.
#' @param pprobs A vector of posterior probabilities for each gene estimated
#' from the
#' intact function. Gene order should match d_vec.
#' @param SE_type A method to compute standard errors of the gene set enrichment
#'  estimates.
#' Possible methods are "profile_likelihood," "bootstrap," and "NDS". NDS
#' performs numerical
#'  differentiation of the Fisher score vector.
#' @param boot_rep Number of bootstrap samples, if boostrap standard errors are
#' specified
#' for SE_type.
#' @return A data frame with the alpha1 estimate, standard error, z-score,
#' p-value,
#' (1-sig_lev)\% CI limits, and convergence indicator.
#' @export
#' @examples
#' data(simdat)
#' enrich_res(d_vec = sample(c(0,1),1197,replace=TRUE),
#' pprobs = intact(GLCP_vec=simdat$GLCP, prior_fun=linear,
#' z_vec = simdat$TWAS_z, t = 0.05),
#'  sig_lev = 0.05)
#' @importFrom numDeriv hessian
#' @importFrom bdsmatrix gchol
#' @importFrom stats pnorm qnorm





enrich_res <- function(sig_lev, pprobs, d_vec, SE_type = "NDS", boot_rep = NULL){

  #compute a1 and a0 MLEs using the SQUAREM package

  alpha <- em_est(pprobs = pprobs, d_vec = d_vec)

  CONVERGED <- alpha[3]

  alpha <- alpha[seq(1,2)]

  #alpha0 MLE

  a0 <- alpha[1]

  #alpha1 MLE

  a1 <- alpha[2]

  #compute the specified standard error type

  if (SE_type == "profile_likelihood"){

    #compute a0 standard error

    SE_a0 <- sqrt(a0^2 / (2*(logistic_loglik(alpha = c(a0,a1),
                                             pprobs = pprobs,
                                             d_vec = d_vec) -
                               logistic_loglik(alpha = c(0,a1),
                                               pprobs = pprobs,
                                               d_vec = d_vec))))


    #compute a1 standard error

    SE_a1 <- sqrt(a1^2 / (2*(logistic_loglik(alpha = c(a0,a1),
                                             pprobs = pprobs,
                                             d_vec = d_vec) -
                               logistic_loglik(alpha = c(a0,0),
                                               pprobs = pprobs,
                                               d_vec = d_vec))))

  }

  if (SE_type == "bootstrap"){

    SE_vec <- enrich_bootstrap_se(pprobs = pprobs, d_vec = d_vec, reps = boot_rep)

    SE_a0 <- SE_vec[1]

    SE_a1 <- SE_vec[2]

  }

  if (SE_type == "NDS"){

    if (NA %in% alpha){

      SE_a0 <- SE_a1 <- NA

    }else{

      #Numerical differentiation of log likelihood function to get Hessian evaluated at MLEs

      hess <- hessian(func = logistic_loglik,
                                x = alpha,
                                method = "Richardson",
                                d_vec = d_vec,
                                pprobs = pprobs)

      #Generalized Cholesky decomposition

      cholesk <- gchol(-1*hess)

      cov_mat <- solve(cholesk)

      SE_a0 <- sqrt(cov_mat[1,1])

      SE_a1 <- sqrt(cov_mat[2,2])

    }
  }
  if (!(SE_type %in% c("profile_likelihood","bootstrap","NDS"))){

    stop("Standard error type must be profile_likelihood, bootstrap, or NDS")

  }

  #assume that likelihood \hat\alpha1|\alpha_1 ~ N(\alpha_1, se^2(\hat alpha_1))
  #and
  #alpha1 ~ N(0,1)
  # so alpha1|\hat\alpha1 ~ N(1/(1+se^2)alpha1, se^2/(se^2 + 1))

  posterior_mean_a1 <- a1*1/(1+SE_a1^2) #shrink a1 estimate

  posterior_var_a1 <- SE_a1^2/(SE_a1^2 + 1)

  posterior_se_a1 <- sqrt(posterior_var_a1) #posterior se for a1

  #If Hessian fails to converge, use a1 prior as the posterior.

  if (posterior_se_a1 %in% c(0,NA)){

    posterior_se_a1 <- 1

    posterior_mean_a1 <- 0

  }

  #If NDS approach fails, use the prior as the posterior.


  out_df <- data.frame(Estimate = c(posterior_mean_a1),
                       SE = c(posterior_se_a1),
                       z = posterior_mean_a1/posterior_se_a1,
                       pval = 2*pnorm(-abs(posterior_mean_a1/posterior_se_a1)),
                       CI_Leftlim = c(posterior_mean_a1 - qnorm(sig_lev/2,
                                                                lower.tail = FALSE)*posterior_se_a1),
                       CI_Rightlim = c(posterior_mean_a1 + qnorm(sig_lev/2,
                                                                 lower.tail = FALSE)*posterior_se_a1),
                       CONVERGED = CONVERGED)

  return(out_df)

}

