% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voomlimmaFit.R
\name{voomlimmaFit}
\alias{voomlimmaFit}
\title{Analysis of RFI RNA-seq data Using voom}
\usage{
voomlimmaFit(counts, design, Effect)
}
\arguments{
\item{counts}{a  matrix of count data.}

\item{design}{a design matrix.}

\item{Effect}{the effect used to simulate data, either line2, or
time. This  effect is considered as the main factor of interest where the
status of DE and EE genes was specified.}
}
\value{
a list of 4 components \item{fit}{output of voom-limma fit.}
  \item{pv}{a vector of p-values of the test for significant of
  \code{Effect}.} \item{qv}{a vector of q-values corresponding to the
  \code{pv} above.}
}
\description{
This function analyzes RFI RNA-seq data and simulated datasets using
\code{\link[limma]{voom}}, which uses precision weights and linear model
pipeline for the analysis of log-transformed RNA-seq data.
}
\examples{
data(dat)
data(design)
counts <- dat[1:50,]
design <- design
Effect <- "line2"
voomlimmaout <- rmRNAseq:::voomlimmaFit(counts, design, Effect)
names(voomlimmaout)
}
\references{
1. Gordon K. Smyth, Matthew Ritchie, Natalie Thorne,James
  Wettenhall, Wei Shi and Yifang Hu. limma: Linear Models for Microarray and
  RNA-Seq Data. User's Guide.
  \url{https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf}

2. Gordon K. Smyth. Linear models and empirical bayes methods for
  assessing differential expression in microarray experiments. Stat Appl
  Genet Mol Biol. 2004;3:Article3. Epub 2004 Feb 12.
}
