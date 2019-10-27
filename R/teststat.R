#' Calculating F-Type Statistics To Test a General Linear Hypothesis
#'
#' This function is to calculate F-type test statistics for a general linear
#' hypothesis for each of G genes.
#'
#' @param C.matrix is a list of matrix Ci in testing H0:  Ci*beta = 0.
#' @param beta0 vector of the hypothesized value of beta, usually, beta0 is a 0
#'   vector. The default option \code{beta0 = NULL} means that \code{beta0} is a
#'   vector of 0.
#' @param regression.output this is a data.frame containing the output of
#'   \code{\link{glsCAR1}} function for all G genes.
#' @param ncores number of cores for embarrassingly parallel procedure. Default
#'   value of \code{ncores} is 1.
#' @return A matrix of dimension G X length(C.matrix) of F-similar test
#'   statistics
#' @examples
#' \donttest{
#' data(design)
#' beta0 <- NULL
#' regression.output <- res$ori.res$newlm[1:50,]
#' ncores <- 1
#' C.matrix <- list()
#' C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
#' C.matrix[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
#' names(C.matrix) <- c("line2","time")
#' teststatout <- rmRNAseq:::teststat(C.matrix, beta0, regression.output, ncores)
#' head(teststatout)
#' }

teststat <- function(C.matrix, beta0=NULL, regression.output, ncores = 1) {
  beta.coef <- data.matrix(regression.output[grep(pattern = "fixed.", x = names(regression.output))])
  lower.tri.vector <- data.matrix(regression.output[grep(pattern = "varbeta", x = names(regression.output))])
  if (is.null(beta0))
    beta0 <- array(0, dim = dim(beta.coef))
  s2_shrunken <- regression.output$s2_shrunken
  s2 <- regression.output$s2
  Ftestout <- parallel::mclapply(1:nrow(regression.output), function(i) {
    m_varb <- varbeta(lower.tri.vector[i, ]) * s2_shrunken[i]/s2[i]
    Fstat <- vapply(C.matrix, function(Cm) as.numeric(crossprod(crossprod(Cm, beta.coef[i,] - beta0[i, ]),
                                                      chol2inv(chol(Matrix::nearPD(crossprod(Cm, m_varb) %*% Cm)$mat))) %*%
                                              crossprod(Cm,beta.coef[i, ] - beta0[i, ])), FUN.VALUE = 1)
    c(Fstat)
  }, mc.cores = ncores)
  Ftestout <- do.call("rbind", Ftestout)
  if (length(C.matrix) == 1) {
    Ftestout <- matrix(Ftestout, ncol = 1)
    colnames(Ftestout) <- paste0("Ftest.", names(C.matrix))
  }
 Ftestout
}
