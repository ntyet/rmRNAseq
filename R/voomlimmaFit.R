#' Analysis of RFI RNA-seq data Using voom
#'
#' This function analyzes RFI RNA-seq data and simulated datasets using
#' \code{\link[limma]{voom}}, which uses precision weights and linear model
#' pipeline for the analysis of log-transformed RNA-seq data.
#' @inheritParams edgeRFit
#' @return a list of 4 components \item{fit}{output of voom-limma fit.}
#'   \item{pv}{a vector of p-values of the test for significant of
#'   \code{Effect}.} \item{qv}{a vector of q-values corresponding to the
#'   \code{pv} above.}
#' @importFrom rlang .data
#' @references 1. Gordon K. Smyth, Matthew Ritchie, Natalie Thorne,James
#'   Wettenhall, Wei Shi and Yifang Hu. limma: Linear Models for Microarray and
#'   RNA-Seq Data. User's Guide.
#'   \url{https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf}
#'
#'
#' @references 2. Gordon K. Smyth. Linear models and empirical bayes methods for
#'   assessing differential expression in microarray experiments. Stat Appl
#'   Genet Mol Biol. 2004;3:Article3. Epub 2004 Feb 12.
#' @examples
#' data(dat)
#' data(design)
#' counts <- dat[1:50,]
#' design <- design
#' Effect <- "line2"
#' voomlimmaout <- rmRNAseq:::voomlimmaFit(counts, design, Effect)
#' names(voomlimmaout)

voomlimmaFit <- function(counts, design, Effect) {
  v <- limma::voom(counts = counts, design = design, lib.size = apply(counts, 2, stats::quantile,
                                                               0.75))
  fit <- limma::lmFit(v)
  if (Effect != "time") {
    C.matrix <- eval(parse(text = paste0("limma::makeContrasts(", Effect, ",levels = design)")))
    fit1 <- limma::contrasts.fit(fit, contrast = C.matrix)
    fit1 <- limma::eBayes(fit1)
    tt <- limma::topTable(fit1, adjust.method = "none", sort.by = "none", number = Inf)
    pv <- tt$P.Value
  }



  if (Effect == "time") {
    # C.matrix <- limma::makeContrasts(.data$time2, .data$time6, .data$time24, levels = design) # for some reason it does not work now!!
    C.matrix <- limma::makeContrasts(time2, time6, time24, levels = design)
    fit1 <- limma::contrasts.fit(fit, contrast = C.matrix)
    fit1 <- limma::eBayes(fit1)
    tt <- limma::topTable(fit1, adjust.method = "none", sort.by = "none", number = Inf)
    pv <- tt$P.Value
  }

  pv <- data.matrix(pv)
  qv <- apply(pv, 2, function(x) jabes.q(x))
  colnames(pv) <- colnames(pv) <- Effect
  res <- list(fit = fit, pv = pv, qv = qv)
  res
}

