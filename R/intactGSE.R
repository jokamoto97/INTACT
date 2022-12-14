#' Perform gene set enrichment estimation and inference, given TWAS scan
#' z-scores and colocalization probabilities.
#'
#' @param gene_data A data frame containing gene names and corresponding
#' colocalization
#'  probabilities and TWAS z-scores for each gene. Column names should be
#'  "gene", "GLCP",
#'  and "TWAS_z'. If the user wishes to specify TWAS Bayes factors instead of
#'  z-scores,
#'  use the column name TWAS_BFs. If the user wishes to specify gene-specific
#'  TWAS priors,
#'  use the column name TWAS_priors.
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
#' intactGSE(gene_data = simdat,gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = step,t = 0.45,
#' gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = expit,t = 0.08,D = 0.08,
#' gene_sets = gene_set_list)
#' intactGSE(gene_data = simdat,prior_fun = hybrid,t = 0.08,D = 0.08,
#' gene_sets = gene_set_list)


intactGSE <- function(gene_data,prior_fun = linear,t = NULL,D = NULL,gene_sets,
                      sig_lev=0.05,SE_type="NDS",boot_rep=NULL){

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
