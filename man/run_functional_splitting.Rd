% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_functional_splitting.R
\name{run_functional_splitting}
\alias{run_functional_splitting}
\title{MAIN FUNCTION: functional splitting}
\usage{
run_functional_splitting(
  scm,
  tfs,
  gda,
  norm = T,
  qtarget = 0.9,
  adjMethod = "bonferroni",
  cutoff = 0.05,
  percExp = 10,
  qExp = 0.05
)
}
\arguments{
\item{scm}{single cell matrix (cells in columns / genes in rows)}

\item{tfs}{vector of TFs names / ID depending on the genes in the scm}

\item{gda}{BP GO file used for the enrichment analysis}

\item{norm}{T = perform seurat normalization (default), F = do not perform the normalization}

\item{qtarget}{top X target genes to consider for each identified cliques}

\item{adjMethod}{multiple test correction to perform, either BH (default) or bonferroni}

\item{cutoff}{p-adjusted value cutoff for the functional enrichment}

\item{percExp}{percentage of cells a gene needs to be expressed in to be considered}

\item{qExp}{quantile of gene expression that will be considered as too lowly expressed and so removed}
}
\value{
functional_split object
}
\description{
MAIN FUNCTION: functional splitting
}
\author{
Celine Barlier
}
