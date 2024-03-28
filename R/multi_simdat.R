#' Simulated TWAS, PWAS, and pairwise colocalization summary data.
#'
#' A data set containing pairwise fastENLOC GLCPs, TWAS and PWAS z-scores, and
#' multivariate Wald chi-square statistics for 1197 simulated genes.
#' @format A data frame with 1197 rows and 6 variables:
#' \describe{
#'   \item{gene}{gene Ensembl ID}
#'   \item{GLCP_1}{Pairwise colocalization probability for the complex trait and
#'   gene expression levels.}
#'   \item{GLCP_2}{Pairwise colocalization probability for the complex trait and
#'   encoded protein levels.}
#'   \item{z_1}{TWAS z-score.}
#'   \item{z_2}{PWAS z-score.}
#'   \item{chisq}{Multivariate Wald chi-square statistic from the regression
#'   of the complex trait on the predicted expression and protein levels of the
#'   target gene.}
#'   ...
#' }
#' @return A data frame with 1197 rows and 6 variables:
#' @examples
#' data(multi_simdat)
"multi_simdat"
