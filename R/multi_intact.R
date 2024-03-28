#' Compute Multi-INTACT prior parameter estimates and gene product relevance
#' probabilities.
#'
#' @param df A data frame containing at least the following six columns:
#' (1)'gene' contains the gene ID,
#' (2)'GLCP_1' contains the numeric pairwise colocalization probability for the
#' complex trait and gene product 1,
#' (3)'GLCP_2' contains the numeric pairwise colocalization probability for the
#' complex trait and gene product 2,
#' (4) 'z_1' A numeric vector of TWAS test statistics for gene product 1,
#' (5) 'z_2' A numeric vector of TWAS test statistics for gene product 2, and
#' (6)'chisq' A numeric vector of multivariate Wald chi-square test statistics.
#' @param chisq_dof A numeric scalar or vector with the degrees of freedom of
#' the chi-square test statistics under the null. Default is 2 for all genes.
#' @param prior_fun A function to transform an aggregated colocalization
#'  probability into a prior.
#' Options are linear, step, expit, and hybrid.
#' @param t A hard threshold for the aggregated colocalization probability.
#'  Values below this number will be
#' shrunk to zero.
#'  This argument is used in the user-specified prior function. Default value
#'  for the step
#'  prior is 0.5. Default value is 0.05 for all other prior functions.
#' @param D A curvature shrinkage parameter. Lower values of D will result in a
#' steeper curve.
#' Default is 0.1. This parameter should only be specified if the user selects
#' the expit or hybrid
#' prior function and does not wish to use the default value.
#' @param xwas_priors An optional vector of user-specified gene-specific priors
#' for the joint regression of the complex trait on the two gene product levels.
#' If no input is supplied, Multi-INTACT computes a scalar prior using the
#' multivariate Wald test data
#' (see the corresponding manuscript for more details).
#' @param xwas_BFs A vector of Bayes factors for the joint regression of the
#' complex trait on the two gene product levels for the gene of interest.
#' This is an alternative
#' option if the user wishes to directly specify Bayes factors instead of
#'  computing them from
#' the chi-square statistics.
#' @param bf_type Method used to compute Bayes factors. Default is 'wakefield',
#' but a formula from Johnson (2005) 'johnson' is also available.
#' @param K A vector of values over which Bayesian model averaging is performed.
#' @param glcp_aggreg Function used to aggregate pairwise coloclization
#' probabilities. Default is taking the max ('max'), but 1-prod(1-GLCP)
#' 'one_minus_prod_one_minus' is also available.
#' @param em_algorithm Should the function run the EM algorithm to compute
#' gene product relevance probabilities (TRUE/FALSE)?
#' @param pi0 Estimate for the parameter pi0.
#' @param pi_init Initialization of prior parameters to be estimated. The first
#' entry should be the qvalue estimate for pi0. The second parameter should
#' be the gene-product-1-only parameter; the third should be the
#' gene-product-2-only parameter, and the fourth should be the
#' gene-product-1+2-parameter. Entries must sum to 1.
#' @param return_model_posteriors Should the function return model posteriors
#' in addition to the gene product relevance probabilities (TRUE/FALSE)?
#' @return A list containing:
#' (1) a data frame with gene probabilities of putative causality (GPPC), model
#' posteriors (posterior_1, posterior_2, and posterior_12) and gene product
#' relevance probabilities (GPRP_1 and GPRP_2);
#' (2) EM algorithm prior parameter estimates;
#' (3) a Boolean indicating convergence of the EM algorithm.
#' @importFrom stats qnorm pchisq
#' @export
#' @examples
#' data(multi_simdat)
#' multi_intact(df = multi_simdat)



multi_intact <- function(df,
                         chisq_dof = 2,
                         prior_fun = linear,
                         t = 0.05, D = NULL,
                         xwas_priors = .pi1_fun(z_vec=qnorm(
                           pchisq(multi_simdat$chisq,df = chisq_dof,
                                  lower.tail = FALSE)/2)),
                         xwas_BFs = NULL,
                         bf_type = "wakefield",
                         K = c(1,2,4,8,16),
                         glcp_aggreg = "max",
                         em_algorithm = TRUE,
                         pi0 = 1 - .pi1_fun(z_vec=qnorm(
                           pchisq(multi_simdat$chisq,df = chisq_dof,
                                  lower.tail = FALSE)/2)),
                         pi_init = c(pi0,rep(1-pi0,3)/3),
                         return_model_posteriors = FALSE){

  if((!("gene" %in% colnames(df)) |
       !("GLCP_1" %in% colnames(df)) |
       !("GLCP_2" %in% colnames(df)) |
       !("z_1" %in% colnames(df)) |
       !("z_2" %in% colnames(df)) |
       !("chisq" %in% colnames(df))) &
     em_algorithm == TRUE &
     is.null(xwas_BFs) == TRUE){

    stop("colnames(df) must contain
         'gene','GLCP_1','GLCP_2','z_1','z_2',and 'chisq'")

  }
  if((!("gene" %in% colnames(df)) |
      !("GLCP_1" %in% colnames(df)) |
      !("GLCP_2" %in% colnames(df)) |
      !("chisq" %in% colnames(df))) &
     em_algorithm == FALSE &
     is.null(xwas_BFs) == TRUE){

    stop("colnames(df) must contain
         'gene','GLCP_1','GLCP_2', and 'chisq'")

  }
  if(length(chisq_dof) > 1 & length(chisq_dof) != nrow(df)){

    stop("The length of chisq_dof must match the number of rows in df.")

  }

  coloc_rst <- df[,c("gene","GLCP_1","GLCP_2")]

  if (glcp_aggreg == "max"){

    coloc_rst$glcp_aggregated <- do.call(pmax, coloc_rst[,2:ncol(coloc_rst)])

  }

  if (glcp_aggreg == "one_minus_prod_one_minus"){

    coloc_rst$glcp_aggregated <- apply(coloc_rst[,2:ncol(coloc_rst)],
                                       MARGIN = 1,
                                       FUN = function(x){
                                         out <- 1-prod(1-x)
                                         return(out)})

  }

  if (!(glcp_aggreg %in% c("max","one_minus_prod_one_minus"))){

    stop("glcp_aggreg must be either max or one_minus_prod_one_minus")

  }


  #Compute prior from GLCPs

  if (length(t) == 0 & length(D) == 0){

    prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,u = xwas_priors)

  }

  if (length(t) == 0 & length(D) > 0){

    prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,u = xwas_priors)

  }

  if (length(t) > 0 & length(D) == 0){

    prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,t = t,u = xwas_priors)

  }

  if (length(t) > 0 & length(D) > 0){

    prior <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,t = t,
                       u = xwas_priors)

  }

  if ("chisq" %in% colnames(df)){

    coloc_rst$chisq_dof <- chisq_dof

    coloc_rst$xwas_pval <- pchisq(df$chisq,
                                  df = coloc_rst$chisq_dof,
                                  lower.tail = FALSE)

    if (bf_type == "johnson"){

      coloc_rst$scan_log10bf <- vapply(df$chisq,
                                       FUN = function(chisq,df){
                                         ifelse(chisq >= df,
                                           log10_bf <- (df/2) *
                                             log10(df/chisq) +
                                             log10(exp((chisq - df)/2)),
                                           log10_bf <- 0)
                                         return(log10_bf)},
                                       FUN.VALUE = numeric(1),
                                       df = coloc_rst$chisq_dof)

    }
    if (bf_type == "wakefield"){

      z_vec <- qnorm(coloc_rst$xwas_pval/2, lower.tail = FALSE)

      bf_grid_log10 <- vapply(X=matrix(K),
                              FUN.VALUE=numeric(length(z_vec)),
                              FUN = function(K,z_vec){
                                out <- 0.5*log10(1/(1+K)) +
                                  (0.5*(z_vec^2) * (K/(1+K)))*log10(exp(1))
                                return(out)
                              },
                              z_vec=z_vec)

      coloc_rst$scan_log10bf <- apply(bf_grid_log10,FUN=function(bf_vec){
        x_star <- max(bf_vec)
        logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
        out <- logsumexp - log10(length(bf_vec))
        return(out)
      },
      MARGIN = 1)

    }
    if (!(bf_type %in% c("johnson","wakefield"))){

      stop("bf_type must be either wakefield or johnson")

    }

    #compute log10(posterior odds)

    log10odds <- coloc_rst$scan_log10bf + log10((prior/(1-prior)))

    #compute GPPCs

    post <- 1/ (1 + 10^(-log10odds))

    coloc_rst$GPPC <- post

    out <- coloc_rst[order(coloc_rst$GPPC,decreasing=TRUE),]

    out <- out[,c("gene","GPPC")]

  }else{ #If XWAS Bayes factors are supplied, compute GPPCs directly.

    #Compute prior from GLCPs

    #compute log10(posterior odds)

    log10odds <- log10(xwas_BFs) + log10((prior/(1-prior)))

    #compute GPPC

    post <- 1/ (1 + 10^(-log10odds))

    coloc_rst$GPPC <- post

    out <- coloc_rst[order(coloc_rst$GPPC,decreasing=TRUE),]

    out <- out[,c("gene","GPPC")]
  }

  #EM algorithm

  if(em_algorithm == FALSE){

    return(out)

  }else{

    if (length(t) == 0 & length(D) == 0){

      fp_coloc <- prior_fun(GLCP = coloc_rst$glcp_aggregated,u = 1)

    }

    if (length(t) == 0 & length(D) > 0){

      fp_coloc <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,u = 1)

    }

    if (length(t) > 0 & length(D) == 0){

      fp_coloc <- prior_fun(GLCP = coloc_rst$glcp_aggregated,t = t,u = 1)

    }

    if (length(t) > 0 & length(D) > 0){

      fp_coloc <- prior_fun(GLCP = coloc_rst$glcp_aggregated,D = D,t = t,
                         u = 1)

    }


    em_rst <- .multi_prior_estimation(df = df,
                                      pi_init = pi_init,
                                      chisq_vec = df$chisq,
                                      chisq_dof = chisq_dof,
                                      z_1 = df$z_1,
                                      z_2 = df$z_2,
                                      fp_coloc = fp_coloc)

    #Merge GPPC results with EM algorithm results

    em_rst[[1]] <- merge(out,em_rst[[1]],by= "gene")

    em_rst[[1]] <- em_rst[[1]][order(em_rst[[1]]$GPPC,decreasing=TRUE),]

    #compute h estimates

    h1 <- em_rst[[2]][2]/(1-em_rst[[2]][1])

    h2 <- em_rst[[2]][3]/(1-em_rst[[2]][1])

    h12 <- em_rst[[2]][4]/(1-em_rst[[2]][1])

    em_rst[[2]] <- c(h1,h2,h12)

    if(return_model_posteriors == TRUE){

    em_rst[[1]] <- em_rst[[1]][,c("gene","GPPC","posterior_0",
                                  "posterior_1","posterior_2","posterior_12",
                                  "GPRP_1","GPRP_2")]
    }else{

      em_rst[[1]] <- em_rst[[1]][,c("gene","GPPC",
                                    "GPRP_1","GPRP_2")]

    }


    names(em_rst) <- c("Muli-INTACT results",
                       "Prior parameter estimates",
                       "EM algorithm converged?")

    return(em_rst)
  }
}

