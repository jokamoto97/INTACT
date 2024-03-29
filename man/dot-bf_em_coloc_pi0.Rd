% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_fixptfn.R
\name{.bf_em_coloc_pi0}
\alias{.bf_em_coloc_pi0}
\title{Multi-INTACT EM algorithm fixed-point function.}
\usage{
.bf_em_coloc_pi0(w, bf, fp_coloc)
}
\arguments{
\item{w}{A vector of prior parameter estimates.}

\item{bf}{A vector of Bayes factors.}

\item{fp_coloc}{A vector of aggregated colocalization evidence. The default
is the max of all GLCPs, truncated at 0.05.}
}
\value{
A vector of updated parameter estimates.
}
\description{
Multi-INTACT EM algorithm fixed-point function.
}
