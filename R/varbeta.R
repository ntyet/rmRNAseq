#' Recovering a Symmetric Matrix from Its Lower Triangular Matrix
#'
#' This function recovers a symmetric matrix from the vector consisting of its
#' lower triangular elements.
#'
#' \code{diag = FALSE} is the case where a correlation matrix whose diagonal
#' elements are 1, therefore there is no need in storing the diagonal elements.
#' @param lower.tri.vector a numerical vector containing the lower triangular
#'   elements of the symmetric matrix.
#' @param diag a logical indicator, \code{TRUE} or \code{FALSE}, indicating if
#'   \code{lower.tri.vector} contains diagonal elements of the symmetric matrix or not.
#' @return The function \code{varbeta} returns the original symmetric matrix. We
#'   only use \code{diag = FALSE} to recover the correlation matrix, because
#'   diagonal elements of a correlation matrix are all 1, so there is no reason
#'   to store these diagonal elements in this case.
#' @author Yet Nguyen \email{ynguyen@odu.edu}
#' @references Yet Nguyen, Dan Nettleton, 2019. rmRNAseq: RNA-seq Analysis for Repeated-measures Data
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' lower.tri.vector <- runif(10)
#' varout <- rmRNAseq:::varbeta(lower.tri.vector, diag=TRUE)
#' varout
#' }
varbeta <- function(lower.tri.vector, diag = TRUE) {
  n <- length(lower.tri.vector)
  if (diag) {
    p <- round((-1 + sqrt(1 + 8 * n))/2)
    X <- diag(p)
    X[lower.tri(X, diag = TRUE)] <- lower.tri.vector
  } else {
    p <- round((1 + sqrt(1 + 8 * n))/2)
    X <- diag(p)
    X[lower.tri(X, diag = FALSE)] <- lower.tri.vector
  }
  X <- X + t(X) - diag(diag(X))
  X
}
