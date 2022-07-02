#' Compute gene set enrichment estimates.
#'
#' @param d_vec A vector of gene set annotations for the genes of interest.
#' Entries should be integer(1) if the gene is annotated and integer(0) otherwise.
#' @param pprobs A vector of posterior probabilities for each gene estimated from the
#' intact function. Gene order should match d_vec.
#' @return Maximum likelihood estimates for alpha0 and alpha1; convergence indicator.
#' @export
#' @examples
#' em_est(d_vec = sample(c(0,1),1197,replace=TRUE), pprobs = intact(GLCP_vec=simdat$GLCP,
#' prior_fun=linear, z_vec = simdat$TWAS_z, t = 0.05))



em_est <- function(pprobs, d_vec){

  alpha_start <- c(log(mean(pprobs)/(1-mean(pprobs))),0)

  CONVERGED <- TRUE

  #Compute MLEs
  square_obj  <- try(SQUAREM::squarem(par = alpha_start,
                                      fixptfn = logistic_em_nopseudo,
                                      objfn = logistic_loglik,
                                      control = list(tol = 1.e-08, minimize=F, maxiter=50),
                                      pprobs=pprobs,
                                      d_vec=d_vec),silent=T)
  if("try-error" %in% class(square_obj)){
    try(square_obj <- SQUAREM::squarem(par = alpha_start,
                                       fixptfn = logistic_em,
                                       objfn = logistic_loglik,
                                       control = list(tol = 1.e-08, minimize=F, maxiter=5),
                                       pprobs=pprobs,
                                       d_vec=d_vec),silent = T)
    CONVERGED <- FALSE
    if("try-error" %in% class(square_obj)){
      return(c(NA,NA,CONVERGED))
    }
    else{
      return(c(square_obj$par,CONVERGED))
    }
  }else{
    return(c(square_obj$par,CONVERGED))
  }
}
