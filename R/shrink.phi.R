#' Shrinkaged Estimates of Error Variance
#'
#' This function implements Smyth's approach (empirical bayes estimate of error
#' variance of genes, limma paper 2004).
#'
#' @param phi.hat a numerical vector of the estimated  of error variance of all
#'   genes.
#' @param den.df denominator degree of freedom associated with the estimated
#'   variances phi.hat (=n sample -rank(design)).
#' @return a list of 3 components \item{phi.shrink}{vector of shrinkaged
#'   estimates of variance.} \item{d0}{estimated prior degree of freedom used
#'   inthe  shrinkage procedure.} \item{phi0}{estimated prior variance used in
#'   the shrinkage procedure.}
#' @examples
#' phi.hat <- rchisq(1000, 1)
#' den.df <- 2
#' shrinkout <- rmRNAseq:::shrink.phi(phi.hat, den.df)
#' hist(shrinkout$phi.shrink)

shrink.phi <- function(phi.hat, den.df) {
  phi.hat[phi.hat <= 0] <- min(phi.hat[phi.hat > 0])
  z <- log(phi.hat)
  z[z == Inf] <- max(z[z != Inf])
  z[z == -Inf] <- min(z[z != -Inf])
  mnz <- mean(z)

  ## solve for d0 and phi0
  d0arg <- stats::var(z) - trigamma((den.df)/2)
  if (d0arg > 0) {
    dif <- function(x, y) abs(trigamma(x) - y)
    inverse.trigamma <- function(y) stats::optimize(dif, interval = c(0, 10000), y = y)$minimum
    d0 <- 2 * inverse.trigamma(d0arg)
    phi0 <- exp(mnz - digamma((den.df)/2) + digamma(d0/2) - log(d0/(den.df)))

    ## compute shrunken phi's
    phi.shrink <- ((den.df) * phi.hat + d0 * phi0)/(den.df + d0)
  } else {
    phi.shrink <- rep(exp(mnz), length(z))
    d0 <- Inf
    phi0 <- exp(mnz)
  }
  return(list(phi.shrink = phi.shrink, d0 = d0, phi0 = phi0))
}
