#' Differential expression analysis based on natural cubic spline regression models for time-course data
#'
#' This function is a modified version of \code{\link[splineTimeR]{splineDiffExprs}} function
#' that allows data input as  log-transformed counts taking into account the voom-weights and 75th quantile as library size.
#' The function compares time dependent behaviour of genes in two different groups.
#' Applying empirical Bayes moderate F-statistic on differences in coefficients of
#' fitted natural cubic spline regression models, differentially expressed in time genes are determined.
#' The function is a wrapper of other R-functions to simplify differential expression analysis of time-course data.
#' @inheritParams splineTimeR::splineDiffExprs
#' @param voom_method logical value TRUE or FALSE. TRUE when using log-transformed counts, voom-weight, and 75th quantile library size
#'
#' @return a list of 4 components
#' \item{fit}{output of \code{\link[limma]{lmFit}} function.}
#' \item{pv}{a vector of p-values of the test for line effect}
#' \item{qv}{a vector of q-values corresponding to the \code{pv} above.}
#' \item{diffExprs}{A data.frame with rows defining names/IDs of differentially expressed genes}
#' @importFrom rlang .data
#' @import splineTimeR
#' @importFrom Biobase fData
#' @importFrom Biobase pData
#' @importFrom Biobase ExpressionSet
#' @examples
#' library(rmRNAseq)
#' data(design)
#' data(covset)
#' data("resSymm")
#' EE <- 40
#' DE <- 10
#' n_idx <- sample(nrow(resSymm$ori.res$v), size = EE + DE)
#' v <- resSymm$ori.res$v[n_idx,]
#' newlm <- resSymm$ori.res$newlm[n_idx,]
#' BetaMat <- data.matrix(newlm[grep("fixed.", names(newlm))])
#' BetaMat[1:EE, 2] <- 0
#' BetaMat[(EE+1):(EE+DE), 2] <- BetaMat[(EE+1):(EE+DE), 2]*4
#' Sigma2Vec <- newlm$s2_shrunken
#' RhoVec <- data.matrix(newlm[grep("rho.", names(newlm))])
#' WeightMat <- v$weights
#' lib.size <- v$targets$lib.size
#' nrep <- 1
#' Subject <- covset$ear
#' Time <- covset$time
#' counts <- rmRNAseq:::sc_Symm(BetaMat, Sigma2Vec, RhoVec, WeightMat,
#' lib.size, design, Subject, Time,nrep)
#' inData <- counts
#' colnames(inData) <- paste(covset$line, Subject, covset$timef,sep="_")
#' design <- rmRNAseq::design
#' design1 <- data.frame(row.names=colnames(inData),
#' "SampleName"=colnames(inData),
#' "Time"=covset$time,"Treatment"=covset$line)
#' phenoData <- new("AnnotatedDataFrame",data=design1)
#' data <- Biobase::ExpressionSet(assayData=as.matrix(inData),phenoData=phenoData)
#' diffExprs <- rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
#'                             reference = "L",intercept = TRUE, voom_method = FALSE)
#'rmRNAseq:::pauc_out(diffExprs$pv, EE , DE)
#'diffExprs2 <- rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
#'reference = "L",intercept = TRUE, voom_method = TRUE)
#'rmRNAseq:::pauc_out(diffExprs2$pv, EE , DE)

my_splineDiffExprs <- function (eSetObject, df, cutoff.adj.pVal = 1, reference, intercept = TRUE, voom_method = FALSE)
{
  if (!is(eSetObject, "ExpressionSet"))
    stop("eSetObject must be of class ExpressionSet")
  if (!(("SampleName" %in% names(Biobase::pData(eSetObject))) & ("Time" %in%
                                                        names(Biobase::pData(eSetObject))) & ("Treatment" %in% names(Biobase::pData(eSetObject)))))
    stop("eSetObject has to include SampleName, Time and Treatment columns in phenotypic data")
  if (!(is(cutoff.adj.pVal, "numeric") & (cutoff.adj.pVal >=
                                          0) & (cutoff.adj.pVal <= 1)))
    stop("cutoff.adj.pVal must be numeric between 0 and 1")
  if (!(is(df, "numeric") & (df%%1 == 0) & (df > 0)))
    stop("df must be integer > 0")
  if (!(reference %in% levels(factor(Biobase::pData(eSetObject)$Treatment))))
    stop("define valid reference group")
  if (!is(intercept, "logical"))
    stop("intercept must be boolean data type")
  b_ <- splines::ns(Biobase::pData(eSetObject)$Time, df = df)
  d_ <- factor(Biobase::pData(eSetObject)$Treatment, levels = c(reference,
                                                       setdiff(levels(factor(Biobase::pData(eSetObject)$Treatment)),
                                                               reference)))

  design <- model.matrix(~d_ * b_)

  # add voom method

  if (voom_method){
    eSetObject_voom <- limma::voom(counts = exprs(eSetObject), design = design,
                                   lib.size = apply(exprs(eSetObject), 2, stats::quantile,
                                                    0.75))
    fit <- limma::lmFit(eSetObject_voom, design)
  }else fit <- limma::lmFit(eSetObject, design)

  fit_eBayes <- limma::eBayes(fit)
  if (ncol(Biobase::fData(eSetObject)) == 0) {
    Biobase::fData(eSetObject) <- data.frame(rownames(exprs(eSetObject)))
    colnames(Biobase::fData(eSetObject)) <- "row_IDs"
  }
  if (intercept) {
    fit_coeff_ref <- fit$coefficient[, c(1, 3:(df + 2))]
    colnames(fit_coeff_ref)[1] <- "b_0"
    Biobase::fData(eSetObject) <- cbind(Biobase::fData(eSetObject), fit_coeff_ref)
    topTable_output <- limma::topTable(fit_eBayes, coef = c(2, (df +
                                                           3):(2 * df + 2)), sort.by = "none", number = Inf,
                                genelist = Biobase::fData(eSetObject))
    colnames(topTable_output)[(ncol(topTable_output) - 4 -
                                 df):(ncol(topTable_output) - 4)] <- paste("d_", 0:df,
                                                                           sep = "")
  }else {
    fit_coeff_ref <- fit$coefficient[, c(1, 3:(df + 2), 2)]
    colnames(fit_coeff_ref)[1] <- "b_0"
    Biobase::fData(eSetObject) <- cbind(Biobase::fData(eSetObject), fit_coeff_ref)
    topTable_output <- limma::topTable(fit_eBayes, coef = c((df +
                                                        3):(2 * df + 2)), sort.by = "none", number = Inf,
                                genelist = fData(eSetObject))
    colnames(topTable_output)[(ncol(topTable_output) - 4 -
                                 df):(ncol(topTable_output) - 4)] <- paste("d_", 0:df,
                                                                           sep = "")
  }
  diffExprs <- topTable_output[topTable_output$adj.P.Val <=
                                 cutoff.adj.pVal, ]

  # diffExprs <- diffExprs[order(diffExprs$adj.P.Val), ] # no ordering the genes
  cat("-------------------------------------------------",
      "\n")
  cat(paste("Differential analysis done for df = ", df, " and adj.P.Val <= ",
            cutoff.adj.pVal, sep = ""), "\n")
  cat("Number of differentially expressed genes: ", nrow(diffExprs),
      "\n")
  pv = topTable_output$P.Value
  qv <- jabes.q(pv)
  out <- list(fit = fit, pv = pv, qv = qv, diffExprs = diffExprs)
  return(out)
}
