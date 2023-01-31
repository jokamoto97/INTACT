#' Perform gene set enrichment estimation and inference, given TWAS scan
#' z-scores and colocalization probabilities.
#'
#' @param gene_data A data frame containing gene names and corresponding
#' colocalization
#'  probabilities and TWAS z-scores for each gene. Column names should be
#'  "gene", "GLCP",
#'  and "TWAS_z'. If the user wishes to specify TWAS Bayes factors instead of
#'  z-scores,
#'  use the column name "TWAS_BFs". If the user wishes to specify gene-specific
#'  TWAS priors,
#'  use the column name "TWAS_priors".
#' @param prior_fun A function to transform a colocalization probability into a
#' prior.
#' Options are linear, step, expit, and hybrid.
#' @param gene_sets A named list of gene sets for which enrichment is to be
#' estimated.
#' List items should be character vectors of gene IDs. Gene ID format should
#' match the
#'  gene column in gene_data.
#' @param t A hard threshold for the GLCP. Values below this number will be
#' shrunk to zero.
#'  This argument is used in the user-specified prior function. Default value
#'   for the step
#'  prior is 0.5. Default value is 0.05 for all other prior functions.
#' @param D A curvature shrinkage parameter. Lower values of D will result in a
#'  steeper curve.
#' Default is 0.1. This parameter should only be specified if the user selects
#' the expit or hybrid
#' prior function and does not wish to use the default value.
#' @param sig_lev A significance threshold for gene set enrichment hypothesis
#' testing.
#' @param SE_type A method to compute standard errors of the gene set enrichment
#'  estimates.
#'  Possible methods are "profile_likelihood" and "bootstrap."
#' @param boot_rep Number of bootstrap samples.
#' @return A data frame with the alpha1 estimate, standard error, z-score,
#' p-value,
#' (1-sig_lev)\% CI limits, and convergence indicator for each gene set in
#' gene_sets.
#' @export
#' @examples
#' data(simdat)
#' data(gene_set_list)
#' intactGSE(gene_data = simdat,gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = step,t = 0.45,
#' gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = expit,t = 0.08,D = 0.08,
#' gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = hybrid,t = 0.08,D = 0.08,
#' gene_sets = gene_set_list)


intactGSE <- function(gene_data,prior_fun = linear,t = NULL,D = NULL,gene_sets,
                      sig_lev=0.05,SE_type="NDS",boot_rep=NULL){

  if (!is.list(gene_sets)){

    stop("gene_sets must be provided as a list object.")

  }


  if (is.data.frame(gene_data) == FALSE){

    stop("gene_data must be provided a data frame object.")

  }

  if (!("gene" %in% colnames(gene_data))){

    stop("One of the columns of gene_data must be 'gene'.")

  }

  if (!("TWAS_z" %in% colnames(gene_data))
      & !("TWAS_BFs" %in% colnames(gene_data))){

    stop("One column of gene_data must be either 'TWAS_z' or 'TWAS_BFs'
         (but both should not exist together in gene_data).")

  }

  if (!("GLCP" %in% colnames(gene_data))){

    stop("One of the columns of gene_data but be 'GLCP'")

  }

  if (sum(!is.na(as.numeric(gene_data$GLCP))) < length(gene_data$GLCP)){

    stop("The GLCP column must be numeric")

  }

  if (sum(gene_data$GLCP < 0 | gene_data$GLCP > 1) != 0){

    stop("The GLCP column must be probabilities between 0 and 1.")

  }

  if (sum(!is.na(as.numeric(gene_data$TWAS_z))) < length(gene_data$TWAS_z)){

    stop("The TWAS_z column must be numeric.")

  }

  if (sum(!is.na(as.numeric(gene_data$TWAS_BFs))) < length(gene_data$TWAS_BFs)){

    stop("The TWAS_BFs column must be numeric.")

  }

  if (sum(!is.na(as.numeric(gene_data$TWAS_priors))) <
      length(gene_data$TWAS_priors)){

    stop("The TWAS_priors column must be numeric.")

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

  if (is.numeric(sig_lev) == FALSE | sig_lev > 1 | sig_lev < 0){

    stop("sig_level must be a number between 0 and 1.")

  }

  if (!(SE_type %in% c("profile_likelihood","bootstrap","NDS"))){

    stop("Standard error type must be profile_likelihood, bootstrap, or NDS")

  }

  if (length(boot_rep) > 0 & SE_type != "bootstrap"){

    stop("Number of bootstrap samples should only be specified if the SE_type
         is specified as bootstrap.")

  }

  if(length(boot_rep) > 0){
     if(boot_rep %% 1 != 0 | boot_rep < 0){

       stop("boot_rep must be a positive integer.")

      }
  }


  rst <- data.frame()

  for (i in seq(1,length(gene_sets))){

    gene_data$d_vec <- 0

    gene_data$d_vec[gene_data$gene %in% gene_sets[[i]]] <- 1

    if ("TWAS_z" %in% colnames(gene_data) & !("TWAS_priors" %in% colnames(gene_data))){

      if (length(t) == 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z)

      }
      if (length(t) > 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            t = t)

      }
      if (length(t) == 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            D = D)

      }
      if (length(t) > 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            t = t,
                            D = D)

      }

    }
    if ("TWAS_z" %in% colnames(gene_data) & ("TWAS_priors" %in% colnames(gene_data))){

      if (length(t) == 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) > 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            t = t,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) == 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            D = D,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) > 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            z_vec = gene_data$TWAS_z,
                            t = t,
                            D = D,
                            twas_priors = gene_data$TWAS_priors)
      }
    }
    if ("TWAS_BFs" %in% colnames(gene_data)){

      if (length(t) == 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            twas_BFs = gene_data$TWAS_BFs,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) > 0 & length(D) == 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            twas_BFs = gene_data$TWAS_BFs,
                            t = t,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) == 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            twas_BFs = gene_data$TWAS_BFs,
                            D = D,
                            twas_priors = gene_data$TWAS_priors)
      }
      if (length(t) > 0 & length(D) > 0){

        posterior <- intact(GLCP_vec = gene_data$GLCP,
                            prior_fun = prior_fun,
                            twas_BFs = gene_data$TWAS_BFs,
                            t = t,
                            D = D,
                            twas_priors = gene_data$TWAS_priors)
      }

    }



    out <- enrich_res(sig_lev = sig_lev,
                      pprobs = posterior,
                      d_vec = gene_data$d_vec,
                      SE_type = SE_type,
                      boot_rep = boot_rep)

    out$Gene_Set <- names(gene_sets)[i]

    out <- out[,c(8,seq(1,7))]

    rst <- rbind.data.frame(rst,out)
  }

  return(rst)
}
