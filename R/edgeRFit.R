#' Analysis of the RFI RNA-seq data Using edgeR
#'
#' This function analyzes the RFI RNA-seq data and simulated datasets using
#' \code{\link[edgeR]{glmQLFTest}}, the Quasi-likelihood F-test in \code{edgeR}
#' package.
#'
#' @param counts a  matrix of count data.
#' @param design a design matrix.
#' @param Effect the effect used to simulate data, either line2, or
#'   time. This  effect is considered as the main factor of interest where the
#'   status of DE and EE genes was specified.
#' @return a list of 4 components
#' \item{fit}{output of \code{\link[edgeR]{glmQLFit}} function.}
#' \item{pv}{a vector of p-values of the test for significant of \code{Effect}.}
#' \item{qv}{a vector of q-values corresponding to the \code{pv} above.}
#' @importFrom rlang .data
#' @examples
#' data(dat)
#' data(design)
#' counts <- dat[1:100,]
#' design <- design
#' Effect <- "line2"
#' edgeRout <- rmRNAseq:::edgeRFit(counts, design, Effect)
#' names(edgeRout)
edgeRFit <- function(counts, design, Effect) {
  y <- edgeR::DGEList(counts = counts, lib.size = apply(counts, 2, stats::quantile, 0.75))
  y <- edgeR::estimateDisp(y, design, robust = T)
  fit <- edgeR::glmQLFit(y, design, robust = T)

  if (!Effect %in% c("int", "time")) {
    C.matrix <- eval(parse(text = paste0("limma::makeContrasts(", Effect, ",levels = design)")))
    qlf <- edgeR::glmQLFTest(fit, contrast = C.matrix, poisson.bound = F)
    pv <- qlf$table$PValue
  }

  if (Effect == "time") {
     C.matrix <- limma::makeContrasts(.data$time2, .data$time6, .data$time24, levels = design)
    qlf <- edgeR::glmQLFTest(fit, contrast = C.matrix, poisson.bound = F)
    pv <- qlf$table$PValue
  }
  pv <- data.matrix(pv)
  qv <- apply(pv, 2, function(x) jabes.q(x))
  colnames(pv) <- colnames(pv) <- Effect
  res <- list(fit = fit, pv = pv, qv = qv)
  res
}
