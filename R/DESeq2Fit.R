#' Analysis of LPS RFI RNA-seq data Using DESeq2
#'
#' This function analyzes the LPS RFI RNA-seq data and simulated datasets using
#' \code{\link[DESeq2]{DESeq}}, the Quasi-likelihood F-test in \code{DESeq2}
#' package.
#'
#' @param counts a  matrix of count data.
#' @param design a design matrix.
#' @param Effect the effect used to simulate data, either line2,  or
#'   time. This  effect is considered as the main factor of interest where the
#'   status of DE and EE genes was specified.
#'  @param covset a data frame contain covariate set.
#' @return a list of 4 components
#' \item{fit}{output of \code{\link[DESeq2]{DESeq}} function.}
#' \item{pv}{a vector of p-values of the test for significant of \code{Effect}.}
#' \item{qv}{a vector of q-values corresponding to the \code{pv} above.}
#' @import DESeq2
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#' data(dat)
#' data(design)
#' counts <- dat[1:100,]
#' design <- design
#' Effect <- "line2"
#' DESeq2Fitout <- rmRNAseq:::DESeq2Fit(counts, design, Effect)
#' names(DESeq2Fitout)
#' }
DESeq2Fit <- function(counts,  design, Effect){
  # register(MulticoreParam(detectCores()))
  counts <- data.matrix(counts)
  covset <- data.frame(design)
  rownames(covset) <- colnames(counts)
  repeat{
    dds <- tryCatch(DESeqDataSetFromMatrix(countData = counts, colData = covset, design = ~line2 + time2 + time6 + time24 + linetime2 +linetime6 +linetime24),
                    error = function(e)NULL)
    if(!is.null(dds))break
    counts[which.max(counts)] <- max(counts[-which.max(counts)]) # for some reason if a gene has an extreme big value, DESeqDataSetFromMatrix return NA for the count dueto external function SummarizedExperiment
  }

  full <- design
  if(!Effect %in%c("int", "time")){
    reduced <- model.matrix(stats::as.formula(paste0("~", paste(setdiff(colnames(full)[-1], Effect), collapse = "+"))), data = covset)
    dds <- DESeq(object = dds, test = "LRT", full = full, reduced = reduced, minReplicatesForReplace =Inf)
    out <- results(dds, test = "LRT", cooksCutoff = Inf)
    pv <- out$pvalue
  }

  # int effect ------------
  if(Effect =="int"){
    reduced <- model.matrix(~line2 + time2 + time6 + time24 , data = covset)
    dds <- DESeq(object = dds, test = "LRT", full = full, reduced = reduced, minReplicatesForReplace =Inf)
    out <- results(dds, test = "LRT", cooksCutoff = Inf)
    pv <- out$pvalue
  }


  # time main effect-------------
  if(Effect =="time"){
    reduced <- model.matrix(~line2 + linetime2 + linetime6 + linetime24 , data = covset)
    dds <- DESeq(object = dds, test = "LRT", full = full, reduced = reduced, minReplicatesForReplace =Inf)
    out <- results(dds, test = "LRT", cooksCutoff = Inf)
    pv <- out$pvalue
  }
  pv <- data.matrix(pv)
  qv <- apply(pv, 2, function(x)jabes.q(x))
  colnames(pv) <- colnames(qv) <- Effect

  res <- list(pv = pv, qv = qv)
  res
}

