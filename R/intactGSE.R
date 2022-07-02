#' Perform gene set enrichment estimation and inference, given TWAS scan z-scores and colocalization probabilities.
#'
#' @param gene_data A dataframe containing gene names and corresponding colocalization
#'  probabilities and TWAS z-scores for each gene. Column names should be "gene", "GLCP",
#'  and "TWAS_z'
#' @param prior_fun A function to transform a colocalization probability into a prior.
#' Options are linear, step, expit, and hybrid.
#' @param gene_sets A named list of gene sets for which enrichment is to be estimated.
#' List items should be character vectors of gene IDs. Gene ID format should match the
#'  gene column in gene_data.
#' @param t A hard threshold for the GLCP. Values below this number will be shrunk to zero.
#' This argument is used in the user-specified prior function.
#' @param sig_lev A significance threshold for gene set enrichment hypothesis testing.
#' @param SE_type A method to compute standard errors of the gene set enrichment estimates.
#'  Possible methods are "profile_likelihood" and "bootstrap."
#' @param boot_rep Number of bootstrap samples.
#' @return A data frame with the alpha1 estimate, standard error, z-score, p-value,
#' (1-sig_lev)\% CI limits, and convergence indicator for each gene set in gene_sets.
#' @export
#' @examples
#' intactGSE(gene_data = simdat,prior_fun = linear,gene_sets = gene_set_list)

intactGSE <- function(gene_data,prior_fun,gene_sets,t=0.05,sig_lev=0.05,SE_type="NDS",boot_rep=NULL){

  rst <- data.frame()

  for (i in 1:length(gene_sets)){

    gene_data$d_vec <- 0

    gene_data$d_vec[gene_data$gene %in% gene_sets[[i]]] <- 1

    posterior <- intact(GLCP_vec = gene_data$GLCP,
                        prior_fun = prior_fun,
                        z_vec = gene_data$TWAS_z,
                        t = t)

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
