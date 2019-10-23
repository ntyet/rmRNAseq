#' Estimate Number of True Null Hypotheses Using Histogram-based Method
#'
#' This function estimates the number of true null hypotheses given a vector of
#' p-values using the method of Nettleton et al. (2006) JABES 11, 337-356. The
#' estimate obtained is identical to the estimate obtained by the iterative
#' procedure described by Mosig et al. Genetics 157:1683-1698. The number of
#' p-values falling into B equally sized bins are counted. The count of each bin
#' is compared to the average of all the bin counts associated with the current
#' bins and all bins to its right.  Working from left to right, the first bin
#' that has a count less than or equal to the average is identified. That
#' average is multiplied by the total number of bins to obtain an estimate of
#' m0, the number of tests for which the null hypothesis is true.
#' @param p a numerical vector of p-value
#' @param B number of bin
#' @return The function returns an estimate of m0, the number of tests for which
#'   the null hypothesis is true.
#' @author Dan Nettleton \email{dnett@iastate.edu}
#' @references 1. Dan Nettleton, J. T. Gene Hwang, Rico A. Caldo and Roger P.
#'   Wise. Estimating the Number of True Null Hypotheses from a Histogram of p
#'   Values. Journal of Agricultural, Biological, and Environmental Statistics
#'   Vol. 11, No. 3 (Sep., 2006), pp. 337-356.
#' @references 2. \url{http://www.public.iastate.edu/~dnett/microarray/multtest.txt}.
#' @examples
#' data(res)
#' p <- res$pqvalue$pv$line2
#' m0 <- rmRNAseq:::estimate.m0(p)
#' m0
estimate.m0 = function(p, B = 20) {
    m <- length(p)
    m0 <- m
    bin <- c(-0.1, (1:B)/B)
    bin.counts = rep(0, B)
    for (i in 1:B) {
        bin.counts[i] = sum((p > bin[i]) & (p <= bin[i + 1]))
    }
    tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
    temp <- bin.counts - tail.means
    index <- min((1:B)[temp <= 0])
    m0 <- B * tail.means[index]
    return(m0)
}

#' Q-value Using Histogram-based Method
#'
#' This function computes q-values using the approach of Nettleton et al.
#'(2006) JABES 11, 337-356.
#'@return The function returns a q-value vector of the input p-value vector.
#' @inheritParams estimate.m0
#' @references 1. Dan Nettleton, J. T. Gene Hwang, Rico A. Caldo and Roger P.
#'   Wise. Estimating the Number of True Null Hypotheses from a Histogram of p
#'   Values. Journal of Agricultural, Biological, and Environmental Statistics
#'   Vol. 11, No. 3 (Sep., 2006), pp. 337-356.
#' @author Dan Nettleton \email{dnett@iastate.edu}
#' @references 2. \url{http://www.public.iastate.edu/~dnett/microarray/multtest.txt}.
#' @examples
#' data(res)
#' p <- res$pqvalue$pv$line2
#' q <- rmRNAseq:::jabes.q(p)
#' sum(q <= .05)
jabes.q = function(p, B = 20) {

    m = length(p)
    m0 = estimate.m0(p, B)
    k = 1:m
    ord = order(p)
    p[ord] = (p[ord] * m0)/(1:m)
    qval = p
    qval[ord] = rev(cummin(rev(qval[ord])))
    return(qval)
}

