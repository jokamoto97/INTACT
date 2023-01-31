#' Compute the posterior probability that a gene may be causal, given a gene's
#' TWAS scan z-score (or Bayes factor) and colocalization probability.
#'
#' @param GLCP_vec A vector of colocalization probabilities for the genes of
#' interest
#' @param prior_fun A function to transform a colocalization probability into
#' a prior.
#' Options are linear, step, expit, and hybrid.
#' @param z_vec A vector of TWAS scan z-scores for the genes of interest.
#' The order of
#' genes must match GLCP_vec.
#' @param t A hard threshold for the GLCP. Values below this number will be
#' shrunk to zero.
#'  This argument is used in the user-specified prior function. Default value
#'  for the step
#'  prior is 0.5. Default value is 0.05 for all other prior functions.
#' @param D A curvature shrinkage parameter. Lower values of D will result in a
#' steeper curve.
#' Default is 0.1. This parameter should only be specified if the user selects
#' the expit or hybrid
#' prior function and does not wish to use the default value.
#' @param K A vector of values over which Bayesian model averaging is performed.
#' @param twas_priors An optional vector of user-specified gene-specific TWAS
#' priors.
#' If no input is supplied, INTACT computes a scalar prior using the TWAS data
#' (see the corresponding manuscript for more details).
#' @param twas_BFs A vector of TWAS Bayes factors for the genes of interest.
#' This is an alternative
#' option if the user wishes to directly specify Bayes factors instead of
#'  computing them from
#' TWAS scan z-scores.
#' @return The vector of posteriors.
#' @export
#' @examples
#' data(simdat)
#' intact(GLCP_vec=simdat$GLCP, z_vec = simdat$TWAS_z)
#' intact(GLCP_vec=simdat$GLCP, prior_fun=expit, z_vec = simdat$TWAS_z,
#' t = 0.02,D = 0.09)
#' intact(GLCP_vec=simdat$GLCP, prior_fun=step, z_vec = simdat$TWAS_z,
#' t = 0.49)
#' intact(GLCP_vec=simdat$GLCP, prior_fun=hybrid, z_vec = simdat$TWAS_z,
#' t = 0.49,D = 0.05)






intact <- function(GLCP_vec, prior_fun = linear, z_vec = NULL, t = NULL,
                   D = NULL, K = c(1,2,4,8,16),
                   twas_priors = pi1_fun(z_vec = z_vec,lambda = 0.5),
                   twas_BFs = NULL){

  if (sum(!is.na(as.numeric(GLCP_vec))) < length(GLCP_vec)){

    stop("GLCP_vec must be a numeric vector.")

  }

  if (sum(GLCP_vec < 0 | GLCP_vec > 1) != 0){

    stop("GLCP_vec must be probabilities between 0 and 1.")

  }

  if (length(GLCP_vec) != length(twas_BFs) & length(GLCP_vec) != length(z_vec)){

    stop("GLCP_vec must be the same length as z_vec or twas_BFs")

  }

  if (sum(!is.na(as.numeric(z_vec))) < length(z_vec)){

    stop("z_vec must be a numeric vector")

  }

  if (length(which(t > 1 | t < 0)) != 0){

    stop("t must be a number between 0 and 1.")

  }

  if (length(t) > 0 & is.numeric(t) == FALSE){

    stop("t must be a number between 0 and 1")

  }

  if (length(D) > 0 & is.numeric(D) == FALSE){

    stop("D must be a number.")

  }

  if (!is.function(prior_fun)){

    stop("Prior function must be one of: linear, step, expit, or hybrid.")

  }

  if (sum(!is.na(as.numeric(K))) < length(K)){

    stop("K must be a numeric vector.")

  }

  if (sum(!is.na(as.numeric(twas_priors))) < length(twas_priors)){

    stop("twas_priors must be a numeric vector.")

  }

  if (sum(!is.na(as.numeric(twas_BFs))) < length(twas_BFs)){

    stop("twas_BFs must be a numeric vector.")

  }


  if (length(z_vec) == 0 & length(twas_BFs) == 0){

    stop("TWAS z-scores or Bayes factors must be supplied.")

  }
  if (length(z_vec) != 0 & length(twas_BFs) != 0){

    stop("Choose to use TWAS z-scores OR Bayes factors, not both.")

  }

  #If z-scores are supplied, compute log10(BF) grid from TWAS z score for each
  #K value

  if (length(z_vec) != 0){

    bf_grid_log10 <- vapply(X=matrix(K),FUN.VALUE=numeric(length(z_vec)),
                           FUN = function(K,z_vec){
      out <- 0.5*log10(1/(1+K)) + (0.5*(z_vec^2) * (K/(1+K)))*log10(exp(1))
      return(out)
    },
    z_vec=z_vec)


    #Perform Bayesian model averaging using log sum exp trick

    bf_log10 <- apply(bf_grid_log10,FUN=function(bf_vec){
      x_star <- max(bf_vec)
      logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
      out <- logsumexp - log10(length(bf_vec))
      return(out)
    },
    MARGIN = 1)

    #Compute prior from GLCPs

    if (length(t) == 0 & length(D) == 0){

      prior <- prior_fun(GLCP = GLCP_vec,u = twas_priors)

    }
    if (length(t) == 0 & length(D) > 0){

      prior <- prior_fun(GLCP = GLCP_vec,D = D,u = twas_priors)

    }
    if (length(t) > 0 & length(D) == 0){

      prior <- prior_fun(GLCP = GLCP_vec,t = t,u = twas_priors)

    }
    if (length(t) > 0 & length(D) > 0){

      prior <- prior_fun(GLCP = GLCP_vec,D = D,t = t,u = twas_priors)

    }
    #compute log10(posterior odds)

    log10odds <- bf_log10 + log10((prior/(1-prior)))

    #compute posterior probability

    post <- 1/ (1 + 10^(-log10odds))

    return(post)

  }else{ #If TWAS Bayes factors are supplied, compute posteriors directly.

    #Compute prior from GLCPs

    if (length(t) == 0 & length(D) == 0){

      prior <- prior_fun(GLCP = GLCP_vec,u = twas_priors)

    }
    if (length(t) == 0 & length(D) > 0){

      prior <- prior_fun(GLCP = GLCP_vec,D = D,u = twas_priors)

    }
    if (length(t) > 0 & length(D) == 0){

      prior <- prior_fun(GLCP = GLCP_vec,t = t,u = twas_priors)

    }
    if (length(t) > 0 & length(D) > 0){

      prior <- prior_fun(GLCP = GLCP_vec,D = D,t = t,u = twas_priors)

    }
    #compute log10(posterior odds)

    log10odds <- log10(twas_BFs) + log10((prior/(1-prior)))

    #compute posterior probability

    post <- 1/ (1 + 10^(-log10odds))

    return(post)
  }
}
