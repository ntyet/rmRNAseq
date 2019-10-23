#' Evaluation of  Differential Expression Analysis Methods
#'
#'This function calculates FDR (false discovery rate), pauc (partial area
#'under ROC curve),  when we know the vector of
#'p-values obtained from a particular differential expression analysis method and the
#'true status of each gene. The function requires that the first EE genes are
#'true EE, and the last DE genes are true DE. This requirement can be fulfilled
#'by reorder the rows of gene expression data set.
#'
#'@param p a vector of p-values.
#'@param EE number of EE genes (the first EE genes in p-value vector p).
#'@param DE number of DE genes (the last DE genes in p-value vector p).
#'@return a vector including  V (the
#'  number of false positives), R (the number of declared positives), FDR (the
#'  false discovery rate), pauc (the partial area under ROC curve with respect
#'  to false positive rate fpr less than or equal to a specified level), and
#'  auc.
#'  Here we consider 3 fpr: 0.05, 0.10, 0.20. So the output includes these 15 elements
#'  and the total auc.
#'
#' @examples
#' set.seed(1)
#' EE <- 1000
#' DE <- 500
#' p1 <- runif(EE)
#' p2 <- rbeta(DE, shape1 = .5, shape2 = 1)
#' p <- c(p1, p2)
#' rmRNAseq:::pauc_out(p, EE, DE)


pauc_out <- function(p, EE, DE) {
  lab <- as.factor(c(rep(0, EE), rep(1, DE)))
  roc.out <- AUC::roc(1 - p, lab)

  fpr <- 0.05
  roc.ind <- sum(roc.out$fpr <= fpr)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- AUC::auc(roc.out, min = roc.min)
  if (length(pauc) == 0)
    pauc <- 0
  qv <- jabes.q(p)
  R <- sum(qv <= fpr)
  V <- sum(which(qv <= fpr) <= EE)
  FDR <- V/max(1, R)
  R.05 <- R; V.05 <- V; FDR.05 <- FDR; pauc.05 <- pauc

  fpr <- 0.10
  roc.ind <- sum(roc.out$fpr <= fpr)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- AUC::auc(roc.out, min = roc.min)
  if (length(pauc) == 0)
    pauc <- 0
  qv <- jabes.q(p)
  R <- sum(qv <= fpr)
  V <- sum(which(qv <= fpr) <= EE)
  FDR <- V/max(1, R)
  R.10 <- R; V.10 <- V; FDR.10 <- FDR; pauc.10 <- pauc


  fpr <- 0.20
  roc.ind <- sum(roc.out$fpr <= fpr)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- AUC::auc(roc.out, min = roc.min)
  if (length(pauc) == 0)
    pauc <- 0
  qv <- jabes.q(p)
  R <- sum(qv <= fpr)
  V <- sum(which(qv <= fpr) <= EE)
  FDR <- V/max(1, R)
  R.20 <- R; V.20 <- V; FDR.20 <- FDR; pauc.20 <- pauc


  auc <- AUC::auc(roc.out)
  out <- c(R.05 = R.05, V.05 = V.05, FDR.05 = FDR.05, pauc.05 = pauc.05,
           R.10 = R.10, V.10 = V.10, FDR.10 = FDR.10, pauc.10 = pauc.10,
           R.20 = R.20, V.20 = V.20, FDR.20 = FDR.20, pauc.20 = pauc.20,
           auc =auc)
  return(out)
}
