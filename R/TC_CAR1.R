
#'Fit General Linear Model with \code{corSymm} Correlation Structure for One
#'Gene
#'
#'This function  fits \code{\link[nlme]{gls}} model with REML estimation method,
#'\code{\link[nlme]{corSymm}} unstructured correlation for one gene in a RNA-seq
#'repeated measures data, where data is the log-transformed counts obtained from
#'\code{\link[limma]{voom}}.
#'
#'@param d a data frame containing several columns. The first 4 columns are
#'  \code{y:} a vector of log-counts (obtained by \code{\link[limma]{voom}}),
#'  \code{Subject:} a vector of subject/experimental units where repeated
#'  measures are obtained (can be either numeric or factor),  \code{Time:} a
#'  vector of time points (continuous, since we fit
#'  \code{\link[nlme]{corSymm}}), \code{w:} weights to put in gls model, this is
#'  the inverse of weights obtained by \code{\link[limma]{voom}} The other
#'  columns are exactly the same as design matrix.
#'@return Output is a vector including the following components
#' \item{aic}{AIC of the fitted model.}
#' \item{s2}{estimate of error variance. }
#' \item{rho}{correlation parameter in \code{corSymm} correlation matrix.}
#' \item{fixed}{fixed effects (estimates of regression parameters).}
#' \item{varbeta}{the estimates of variance of fixed effects, just include
#'  lower part and diagonal part of the variance-covariance matrix.}
#' @examples
#'data(res)
#'data(design)
#'data(covset)
#'d <- data.frame(cbind(y = res$ori.res$v$E[1,] ,Subject = covset$ear,
#'Time = as.integer(as.factor(covset$time)), w = 1/res$ori.res$v$weights[1,], design))
#'glsout <- rmRNAseq:::glsSymm(d)
#'glsout

glsSymm <- function(d) {
  d1 <- d
  cnt <- 0
  repeat {
    cnt <- cnt + 1
    # message(cnt, '\n')
    if (cnt > 100) {
      # message(cnt, '\n')
      stop("Too many loops, cnt >100")
    }
    glsout <- tryCatch(nlme::gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)],
                                                                         collapse = "+"))), data = d, weights = nlme::varFixed(~w), correlation = nlme::corSymm(form = ~ Time |
                                                                                                                                                                  Subject), method = "REML", verbose = FALSE, control = nlme::glsControl(maxIter = 1e+07,
                                                                                                                                                                                                                                         msMaxIter = 1e+07, opt = "nlminb", tolerance = 1e-04)), error = function(e) NULL)

    if (is.null(glsout))
      glsout <- tryCatch(nlme::gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)],
                                                                           collapse = "+"))), data = d, weights = nlme::varFixed(~w), correlation = nlme::corSymm(form = ~ Time |
                                                                                                                                                                    Subject), method = "REML", verbose = FALSE, control = nlme::glsControl(maxIter = 1e+07,
                                                                                                                                                                                                                                           msMaxIter = 1e+07, tolerance = 1e-04, opt = "optim")), error = function(e) NULL)

    if (!is.null(glsout)) {
      break
    } else {
      d$y <- d1$y + pmax(-0.1, pmin(stats::rnorm(nrow(d), 0, 0.1), 0.1))
    }
  }
  aic <- stats::AIC(glsout) - 2*(ncol(d)-4) # 4 is the first 4 columns y, Subject, Time, w that are not in design matrix
  ss <- length(unique(d$Time))*(length(unique(d$Time))-1)/2 +1 # number of within-group variance-covariance parameters
  bic <- stats::AIC(glsout) - 2*(ncol(d)-4 + ss) + ss*log(nrow(d) - ncol(d)+4) #
  s2 <- glsout$sigma^2
  fixed <- stats::coef(glsout)
  temp <- stats::vcov(glsout)
  temp <- Matrix::nearPD(temp)$mat
  rho <- unname(stats::coef(glsout$modelStruct$corStruct, unconstrained = FALSE))
  varbetavec <- temp[lower.tri(temp, diag = TRUE)]
  return(c(aic = aic, bic = bic, s2 = s2, rho = rho, fixed = stats::coef(glsout), varbeta = varbetavec))
}


#' General Linear Model Using Voom Output corSymm correlction structure
#'
#' This function run general linear model with \code{\link[nlme]{corSymm}}
#' correlation structure in function \code{\link[nlme]{gls}} for all genes where
#' the input data come from the output of \code{\link[limma]{voom}}.
#' @inheritParams teststat
#' @param v output of \code{\link[limma]{voom}} function.
#' @param Subject a vector of subjects or experimental units.
#' @param Time a vector of time points.
#' @param print.progress logical indicator, T or F, to print the progress.
#'
#' @return a data frame has G rows (= number of genes) containing all outputs from
#'   \code{\link{glsSymm}} function, shrinkage estimates of error variances, and
#'   F-type test statistics calculated by  \code{\link{teststat}} function.
#' @examples
#' data(res)
#' data(covset)
#' v <- res$ori.res$v
#' v$E <- v$E[1:20,]
#' v$weights <- v$weights[1:20,]
#' Subject <- covset$ear
#' Time <- covset$time
#' ncores <- 1
#' C.matrix <- list()
#' C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
#' C.matrix[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
#' names(C.matrix) <- c("line2", "time")
#' beta0 <- NULL
#' print.progress <- FALSE
#' voomglsout <- rmRNAseq:::voomgls_Symm(v, Subject, Time, ncores, C.matrix, beta0, print.progress)
voomgls_Symm <- function(v, Subject, Time, ncores=1, C.matrix = NULL, beta0=NULL, print.progress = FALSE) {
  design <- v$design
  W <- v$weights
  vout <- if(print.progress){
    parallel::mclapply(1:nrow(v), function(i) {
      # print(i)
      d <- data.frame(cbind(y = v$E[i, ], Subject = Subject, Time = as.integer(as.factor(Time)), w = 1/W[i,
                                                                                  ], design))
      message('\n Analyzing gene = ', i, '\n')
      glsSymm(d)
    }, mc.cores = ncores)
  }else{
    parallel::mclapply(1:nrow(v), function(i) {
      d <- data.frame(cbind(y = v$E[i, ], Subject = Subject, Time = as.integer(as.factor(Time)), w = 1/W[i,
                                                                                  ], design))
      glsSymm(d)
    }, mc.cores = ncores)
  }
  vout <- do.call("rbind", vout)
  vout <- data.frame(vout)
  vout$s2_shrunken <- shrink.phi(vout$s2, den.df = nrow(design) - ncol(design))$phi.shrink
  vout$outlier <- 0
  beta.coef <- data.matrix(vout[, grep(pattern = "fixed.", x = names(vout))])
  if (is.null(beta0))beta0 <- array(0, dim = dim(beta.coef))
  if (!is.null(C.matrix)) {
    Ftest <- teststat(C.matrix, beta0 = beta0, regression.output = vout, ncores = ncores)
    vout <- cbind(vout, Ftest = Ftest)
  }
  vout
}


#'Simulating Count Data From The Output of Real Data Analysis (corSymm)
#'
#'This function generates bootstrap samples using parametric bootstrap method.
#'
#'@param BetaMat a matrix of estimates of regression coefficients.
#'@param Sigma2Vec a vector of shrinkage estimates of error variances.
#'@param RhoVec a vector of estimates of correlation.
#'@param WeightMat a matrix of weights of all genes obtaining from voom.
#'@param lib.size library size in voom method, we choose .75 quantile as library
#'  size.
#'@param design a design matrix.
#'@param Subject a vector of subjects/experimental units.
#'@param Time a vector of time points.
#'@param nrep simulation iteration.
#'@return a matrix of count data that has nrow(BetaMat) rows and nrow(design)
#'  columns.

#' @examples
#' \donttest{
#'data(resSymm)
#'v <- resSymm$ori.res$v[1:20,]
#'newlm <- resSymm$ori.res$newlm[1:20,]
#'BetaMat <- data.matrix(newlm[grep("fixed.", names(newlm))])
#'Sigma2Vec <- newlm$s2_shrunken
#'RhoVec <- data.matrix(newlm[grep("rho.", names(newlm))])
#'WeightMat <- v$weights
#'lib.size <- v$targets$lib.size
#'nrep <- 1
#'Subject <- covset$ear
#'Time <- covset$time
#'simcounts <- rmRNAseq:::sc_Symm(BetaMat, Sigma2Vec, RhoVec, WeightMat,
#'lib.size, design, Subject, Time,nrep)
#'dim(simcounts)
#'}

sc_Symm <- function(BetaMat, Sigma2Vec, RhoVec, WeightMat, lib.size, design, Subject, Time, nrep) {
  set.seed(10^9 + nrep)
  Time <- as.integer(as.factor(Time))
  n <- length(Time)
  nSubject <- length(unique(Subject))
  uTime <- unique(Time)
  nTime <- length(uTime)
  # Time <- as.integer(as.factor(Time))
  sim.counts <- vapply(1:nrow(BetaMat), function(i) {
    Xbeta <- design %*% as.numeric(BetaMat[i, ])
    # block.matrix <- Sigma2Vec[i] * corr_matrix(RhoVec[i,], diag = FALSE)
    # block.matrix <- Matrix::nearPD(block.matrix)$mat  # library(Matrix)
    # block.matrix <- matrix(drop(block.matrix@x), nrow = nTime, byrow = TRUE)
    # V <- diag(sqrt(1/WeightMat[i, ])) %*% kronecker(diag(nSubject), block.matrix) %*%
    #   diag(sqrt(1/WeightMat[i, ]))
    corr.matrix <- corr_matrix(RhoVec[i,], diag = FALSE)
    corr.matrix <- Matrix::nearPD(corr.matrix)$mat  # library(Matrix)
    corr.matrix <- matrix(drop(corr.matrix@x), nrow = nTime, byrow = TRUE)
    block.matrix <- matrix(0, nrow = n, ncol = n)
    for(i1 in 1:n){
      for(j1 in 1:n)
        if(Subject[i1] == Subject[j1]) block.matrix[i1,j1] <-corr.matrix[Time[i1], Time[j1]]
    }
    V <- Sigma2Vec[i]*diag(sqrt(1/WeightMat[i, ])) %*%block.matrix%*%diag(sqrt(1/WeightMat[i, ]))
    cnt <- 1
    repeat {
      epsilon <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = V)  # library(MASS)
      out <- unname(round(2^(Xbeta + epsilon) * (lib.size + 1)/10^6)[, , drop = TRUE])
      if (mean(out)>1 & max(out) < 10^7)
        break
    }
    out
  }, FUN.VALUE = rep(1, nrow(design)))
  t(sim.counts)
}



#' myvoom function
#'
#' This function modifies the original \code{\link[limma]{voom}} function
#' to obtain the fitted function f of the lowess fit.
#'
#' @inheritParams limma::voom
#' @importFrom methods is new
#' @importFrom stats approxfun as.formula lowess
#' @importFrom graphics lines title
#' @importFrom Biobase fData pData exprs
#' @return the same output as \code{\link[limma]{voom}} plus the fitted function f
#' of the lowess fit of the \code{link[limma]{voom}} method.
#'
myvoom <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none",
                    span = 0.5, plot = FALSE, save.plot = FALSE, ...) {
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
        0)
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size))
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts)))
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts)))
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size))
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- limma::normalizeBetweenArrays(y, method = normalize.method)
  fit <- limma::lmFit(y, design, ...)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets))
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )",
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  out$f <- f
  out$fit <- fit
  new("EList", out)
}


#'Fit General Linear Model with \code{corCAR1} Correlation Structure for One
#'Gene
#'
#'This function  fits \code{\link[nlme]{gls}} model with REML estimation method,
#'\code{\link[nlme]{corCAR1}} correlation structure for one gene in a RNA-seq
#'repeated measures data, where data is log-transformed output from
#'\code{\link[limma]{voom}}.
#'
#'@param d a data frame containing several columns. The first 4 columns are
#'  \code{y:} a vector of log-counts (obtained by \code{\link[limma]{voom}}),
#'  \code{Subject:} a vector of subject/experimental units where repeated
#'  measures are obtained (can be either numeric or factor),  \code{Time:} a
#'  vector of time points (continuous, since we fit
#'  \code{\link[nlme]{corCAR1}}), \code{w:} weights to put in gls model, this is
#'  the inverse of weights obtained by \code{\link[limma]{voom}} The other
#'  columns are exactly the same as design matrix.
#'@return Output is a vector including the following components
#' \item{aic}{AIC of the fitted model.}
#' \item{s2}{estimate of error variance. }
#' \item{rho}{correlation parameter in \code{corCAR1} correlation matrix.}
#' \item{fixed}{fixed effects (estimates of regression parameters).}
#' \item{varbeta}{the estimates of variance of fixed effects, just include
#'  lower part and diagonal part of the variance-covariance matrix.}
#' @examples
#' \donttest{
#'data(res)
#'data(design)
#'data(covset)
#'d <- data.frame(cbind(y = res$ori.res$v$E[1,] ,Subject = covset$ear,
#'Time = covset$time, w = 1/res$ori.res$v$weights[1,], design))
#'glsout <- rmRNAseq:::glsCAR1(d)
#'glsout
#'}

glsCAR1 <- function(d) {
  d1 <- d
  cnt <- 0
  repeat {
    cnt <- cnt + 1
    # message(cnt, '\n')
    if (cnt > 100) {
      # message(cnt, '\n')
      stop("Problem with glsCAR1, too many loops, cnt >100")
    }
    glsout <- tryCatch(nlme::gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)],
                                                                         collapse = "+"))), data = d, weights = nlme::varFixed(~w), correlation = nlme::corCAR1(form = ~Time |
                                                                                                                                                                  Subject), method = "REML", verbose = FALSE, control = nlme::glsControl(maxIter = 1e+07,
                                                                                                                                                                                                                                         msMaxIter = 1e+07, opt = "nlminb", tolerance = 1e-04)), error = function(e) NULL)

    if (is.null(glsout))
      glsout <- tryCatch(nlme::gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)],
                                                                           collapse = "+"))), data = d, weights = nlme::varFixed(~w), correlation = nlme::corCAR1(form = ~Time |
                                                                                                                                                                    Subject), method = "REML", verbose = FALSE, control = nlme::glsControl(maxIter = 1e+07,
                                                                                                                                                                                                                                           msMaxIter = 1e+07, tolerance = 1e-04, opt = "optim")), error = function(e) NULL)

    if (!is.null(glsout)) {
      break
    } else {
      d$y <- d1$y + pmax(-0.1, pmin(stats::rnorm(nrow(d), 0, 0.1), 0.1))
    }
  }
  aic <- stats::AIC(glsout) - 2*(ncol(d)-4)
  ss <- 1+1
  bic <- stats::AIC(glsout) - 2*(ncol(d)-4 + ss) + ss*log(nrow(d) - ncol(d)+4)
  s2 <- glsout$sigma^2
  rho <- unname(stats::coef(glsout$modelStruct$corStruct, unconstrained = FALSE))
  fixed <- stats::coef(glsout)
  temp <- stats::vcov(glsout)
  temp <- Matrix::nearPD(temp)$mat
  varbeta <- temp[lower.tri(temp, diag = TRUE)]  # save vcov matrix by using only lower triangular including diagonal elements
  return(c(aic = aic, bic = bic, s2 = s2, rho = rho, fixed = stats::coef(glsout), varbeta = varbeta))
}


#' General Linear Model Using Voom Output
#'
#' This function run general linear model with \code{\link[nlme]{corCAR1}}
#' correlation structure in function \code{\link[nlme]{gls}} for all genes where
#' the input data come from the output of \code{\link[limma]{voom}}.
#' @inheritParams teststat
#' @param v output of \code{\link[limma]{voom}} function.
#' @param Subject a vector of subjects or experimental units.
#' @param Time a vector of time points.
#' @param print.progress logical indicator, TRUE or FALSE, to print the progress.
#'
#' @return a data frame has G rows (= number of genes) containing all outputs from
#'   \code{\link{glsCAR1}} function, shrinkage estimates of error variances, and
#'   F-type test statistics calculated by  \code{\link{teststat}} function.
#' @examples
#' \donttest{
#' data(res)
#' data(covset)
#' v <- res$ori.res$v
#' v$E <- v$E[1:50,]
#' v$weights <- v$weights[1:50,]
#' Subject <- covset$ear
#' Time <- covset$time
#' ncores <- 1
#' C.matrix <- list()
#' C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
#' C.matrix[[2]] <- limma::makeContrasts(time2, levels = design)
#' C.matrix[[3]] <- limma::makeContrasts(time6, levels = design)
#' C.matrix[[4]] <- limma::makeContrasts(time24, levels = design)
#' C.matrix[[5]] <- limma::makeContrasts(linetime2, levels = design)
#' C.matrix[[6]] <- limma::makeContrasts(linetime6, levels = design)
#' C.matrix[[7]] <- limma::makeContrasts(linetime24, levels = design)
#' C.matrix[[8]] <- limma::makeContrasts(time2, time6, time24, levels = design)
#' C.matrix[[9]] <- limma::makeContrasts(linetime2,linetime6, linetime24, levels = design)
#' names(C.matrix) <- c("line2", "time2", "time6", "time24",
#'                     "linetime2", "linetime6", "linetime24",
#'                     "time", "int")
#' beta0 <- NULL
#' print.progress <- FALSE
#' voomglsout <- rmRNAseq:::voomgls_CAR1(v, Subject, Time, ncores,
#' C.matrix, beta0, print.progress)
#' }
#' @export
voomgls_CAR1 <- function(v, Subject, Time, ncores=1, C.matrix, beta0=NULL, print.progress = FALSE) {
  design <- v$design
  W <- v$weights
  vout <- if(print.progress){
    parallel::mclapply(1:nrow(v), function(i) {
      # print(i)
      d <- data.frame(cbind(y = v$E[i, ], Subject = Subject, Time = Time, w = 1/W[i,
                                                                                  ], design))
      # message('\n Analyzing gene = ', i, '\n')
      glsCAR1(d)
    }, mc.cores = ncores)
  }else{
    parallel::mclapply(1:nrow(v), function(i) {
      d <- data.frame(cbind(y = v$E[i, ], Subject = Subject, Time = Time, w = 1/W[i,
                                                                                  ], design))
      glsCAR1(d)
    }, mc.cores = ncores)
  }
  vout <- do.call("rbind", vout)
  vout <- data.frame(vout)
  vout$s2_shrunken <- shrink.phi(vout$s2, den.df = nrow(design) - ncol(design))$phi.shrink
  if (!is.null(C.matrix)) {
    Ftest <- teststat(C.matrix, beta0 = beta0, regression.output = vout, ncores = ncores)
    vout <- cbind(vout, Ftest = Ftest)
  }
  vout
}


#' Identify Time Points Mapping to 0 and 1
#'
#' This function is to identify which time points
#' mapped to 0 and 1 based on the estimated correlations
#' of observations between all pairs of time points.
#' The correlation parameters are estimateed by fitting
#' \code{\link{voomgls_Symm}} to the \code{\link[limma]{voom}} transformed data.
#' @inheritParams voomgls_Symm
#' @importFrom utils combn
#' @importFrom stats coef quantile vcov
#' @return a vector of 2 components correspinding to the two time points
#' that are mapping to 0, 1, respectively.
#' @examples
#' \donttest{
#' data(res)
#' data(covset)
#' v <- res$ori.res$v
#' v$E <- v$E[1:2,]
#' v$weights <- v$weights[1:2,]
#' Subject <- covset$ear
#' Time <- covset$time
#' ncores <- 1
#' TimeMinout <- rmRNAseq:::TimeMin(v, Subject, Time, ncores)
#' TimeMinout
#' }
TimeMin <- function(v, Subject, Time, ncores){
  voomglsout_Symm <- voomgls_Symm(v, Subject, Time, ncores, C.matrix = NULL, beta0= NULL, print.progress = FALSE)
  RhoVec <- data.matrix(voomglsout_Symm[grep("rho", names(voomglsout_Symm))])
  uTime <- unique(Time)
  rank_uTime <- rank(uTime)
  # uTimePair <- combn(uTime, 2, simplify = FALSE)
  uTimePair <- combn(1:length(uTime), 2, simplify = FALSE) # corSymm index correlation parameter by consecutive integers from 1 to length(uTime)
  out <- apply(RhoVec, 2, median)
  MinMaxTime <- uTime[rank_uTime%in%uTimePair[[which.min(out)]]]
  # MinMaxTime <- c(min(MinMaxTimev), max(MinMaxTimev))

  outall <- list(MinMaxTime = MinMaxTime, medianout = out, voomglsout_Symm = voomglsout_Symm)
  outall
}


#' Calculate REML Log-Likelihood of glsCAR1 model for each gene
#'
#' This function calculates log-likelihood value of glsCAR1 for
#' each gene using \code{\link[limma]{voom}} data.
#'
#' @inheritParams glsCAR1
#' @return reml log-likelihood value
#' @importFrom nlme gls varFixed corCAR1 glsControl
#' @examples
#' \donttest{
#'data(res)
#'data(design)
#'data(covset)
#'d <- data.frame(cbind(y = res$ori.res$v$E[1,] ,Subject = covset$ear,
#'Time = covset$time, w = 1/res$ori.res$v$weights[1,], design))
#'glsloglikout <- rmRNAseq:::glsCAR1_loglik(d)
#'glsloglikout
#'}


glsCAR1_loglik <- function(d){
  d1 <- d
  cnt <- 0
  repeat{
    cnt <- cnt+1
    # message(cnt, "\n")
    if(cnt >100){
      # message(cnt, "\n")
      stop("Problem with glsCAR1, loop take too long time, cnt >100")
    }
    glsout <- tryCatch(nlme::gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)], collapse = "+"))),
                           data = d,
                           weights = varFixed(~w),
                           correlation = corCAR1( form = ~Time|Subject),
                           method = "REML",verbose = FALSE,
                           control = glsControl(maxIter = 10000000,msMaxIter = 10000000, opt="nlminb",
                                                tolerance = 1e-4
                           )), error = function(e)NULL)

    if (is.null(glsout))glsout <- tryCatch(gls(stats::as.formula(paste0("y~0+", paste0(colnames(d)[-c(1:4)], collapse = "+"))),
                                               data = d,
                                               weights = varFixed(~w),
                                               correlation = corCAR1( form = ~Time|Subject),
                                               method = "REML",verbose = FALSE,
                                               control = glsControl(maxIter = 10000000,msMaxIter = 10000000,
                                                                    tolerance = 1e-4, opt = "optim"
                                               )), error = function(e)NULL)
    # if(!is.null(glsout)&!is.character(glsout$apVar)){
    if(!is.null(glsout)){
      break
    }else{
      d$y <- d1$y + pmax(-.1, pmin(stats::rnorm(nrow(d), 0, 0.1), .1))
      # d$y <- d$y +4
    }
  }
  rmle <- as.numeric(logLik(glsout))
  rmle
}

#' Estimate New Time Points
#'
#' This function estimate new time points to fit the glsCAR1
#' for the \code{\link[limma]{voom}} transformed data. Note that
#' this function is very specific for this dataset, with only
#' 4 time points. If there are more than 4 time points, the
#' method needs to be updated.
#' @inheritParams voomgls_CAR1
#' @param TimeMinOut output from the \code{\link{TimeMin}} function
#' @inheritParams TimeMin
#' @import parallel
#' @importFrom stats constrOptim logLik median
#' @return New time points.
#' @examples
#' \donttest{
#' data(res)
#' data(covset)
#' v <- res$ori.res$v
#' v$E <- v$E[1:2,]
#' v$weights <- v$weights[1:2,]
#' Subject <- covset$ear
#' Time <- covset$time
#' ncores <- 1
#' tmOut <- rmRNAseq:::TimeMin(v, Subject, Time, ncores)
#' TimeMinOut <-tmOut$MinMaxTime
#' NewTimeOut <- rmRNAseq:::NewTimeEst(v, Subject, Time, TimeMinOut, ncores)
#' NewTimeOut
#' }
NewTimeEst <- function(v, Subject, Time, TimeMinOut, ncores){
  TimeOther <- setdiff(unique(Time), TimeMinOut)
  Time1 <- Time
  xmap <- function(x){
    #x0 <- x[1]; x2 <- x[2]
    Time1[Time==TimeMinOut[1]] <- 1
    Time1[Time==TimeMinOut[2]] <- 0
    Time1[Time==TimeOther[1]] <- x[1]
    Time1[Time==TimeOther[2]] <- x[2]
    design <- v$design
    W <- v$weights
    vout <- mclapply(1:nrow(v$E), function(i){
      d <- data.frame(cbind(y = v$E[i, ], Subject =Subject, Time = Time1, w = 1/W[i,], design))
      glsCAR1_loglik(d)
    },mc.cores = ncores)
    vout <- do.call("c", vout)
    loglik <- sum(vout)
    -loglik
  }
  ui <- matrix(c(1, 0,
                 0, -1), byrow = TRUE, ncol = 2)
  ci <- c(0,
          -1)
  o <- constrOptim(theta = c(.3, .51), f = xmap,grad = NULL,  ui = ui, ci = ci)
  out <- o$par
  names(out) <- c("t0", "t2")
  out
}



#'Simulating Count Data From The Output of Real Data Analysis
#'
#'This function generates bootstrap samples using parametric bootstrap method.
#'
#'@param BetaMat a matrix of estimates of regression coefficients.
#'@param Sigma2Vec a vector of shrinkage estimates of error variances.
#'@param RhoVec a vector of estimates of correlation.
#'@param WeightMat a matrix of weights of all genes obtaining from voom.
#'@param lib.size library size in voom method, we choose .75 quantile as library
#'  size.
#'@param design a design matrix.
#'@param Subject a vector of subjects/experimental units.
#'@param Time a vector of time points.
#'@param nrep simulation iteration.
#'@return a matrix of count data that has nrow(BetaMat) rows and nrow(design)
#'  columns.
#' @references Yet Nguyen, Dan Nettleton, 2019. rmRNAseq: RNA-seq Analysis for Repeated-measures Data.
#' @export
#' @examples
#' \donttest{
#'data(res)
#'v <- res$ori.res$v[1:50,]
#'newlm <- res$ori.res$newlm[1:50,]
#'BetaMat <- data.matrix(newlm[grep("fixed.", names(newlm))])
#'Sigma2Vec <- newlm$s2_shrunken
#'RhoVec <- newlm$rho
#'WeightMat <- v$weights
#'lib.size <- v$targets$lib.size
#'nrep <- 1
#'Subject <- covset$ear
#'Time <- covset$time
#'simcounts <- rmRNAseq:::sc_CAR1(BetaMat, Sigma2Vec, RhoVec, WeightMat,
#'lib.size, design, Subject, Time, nrep)
#'dim(simcounts)
#'}

sc_CAR1 <- function(BetaMat, Sigma2Vec, RhoVec, WeightMat, lib.size, design, Subject, Time, nrep) {
  set.seed(10^8 + nrep)
  nSubject <- length(unique(Subject))
  uTime <- unique(Time)
  nTime <- length(uTime)
  corematrix <- abs(outer(uTime, uTime, FUN = "-"))
  sim.counts <- vapply(1:nrow(BetaMat), function(i) {
    # cat (i, '\n')
    Xbeta <- design %*% as.numeric(BetaMat[i, ])
    block.matrix <- Sigma2Vec[i] * RhoVec[i]^corematrix
    block.matrix <- Matrix::nearPD(block.matrix)$mat  # library(Matrix)
    block.matrix <- matrix(drop(block.matrix@x), nrow = nTime, byrow = TRUE)
    V <- diag(sqrt(1/WeightMat[i, ])) %*% kronecker(diag(nSubject), block.matrix) %*%
      diag(sqrt(1/WeightMat[i, ]))
    cnt <- 1
    repeat {
      #message(cnt, '\n')
      epsilon <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = V)  # library(MASS)
      # out <- unname(round(2^(Xbeta + epsilon) * (lib.size + 1)/10^6)[, , drop = TRUE])
      out <- round(pmax(0,as.numeric(2^(Xbeta + epsilon) * (lib.size + 1)/10^6)))
      if (mean(out)>1 & max(out) < 10^7)
        break
      #cnt <- cnt+1
    }
    out
  }, FUN.VALUE = rep(1, nrow(design)))
  t(sim.counts)
}

#' RNA-seq Analysis for Repeated-measures Data
#'
#' This function implements our  parametric bootstrap to analyze
#' repeated measures RNA-seq data.
#'
#'
#' @inheritParams voomgls_CAR1
#' @param counts a matrix of RNA-seq counts.
#' @param design a design matrix.
#' @param Subject a vector of subjects or experimental units.
#' @param Time a vector of time points.
#' @param C.matrix is a list of matrix Ci in testing H0:  Ci*beta = 0.
#' @param Nboot number of bootstrap replicates, default is 100.
#' @param saveboot TRUE or FALSE to save or not save bootstrap output
#' @param circadian TRUE or FALSE to indicate if there is circadian rhythm
#'
#' @return a list of 3 components \item{ori.res}{a list of 2 components
#'   \code{v}: voom's output, \code{newlm}: output from \code{\link{voomgls_CAR1}}.}
#'   \item{boot.res}{a list of \code{Nboot} components, each component is the
#'   output of \code{\link{voomgls_CAR1}} when apply this function to the
#'   corresponding bootstrap sample.} \item{pqvalue}{a list 2 components:
#'   \code{pv}: a matrix of p-values of the tests construted in C.matrix
#'   \code{qv}: matrix of q-values obtaining from using Nettleton 2006 paper
#'   approach, using \code{\link{jabes.q}} function.}
#' @references Yet Nguyen, Dan Nettleton, 2019. rmRNAseq: RNA-seq Analysis for Repeated-measures Data.
#' @examples
#' # This example shows how to implement the method using LPS RFI data.
#' data(dat)
#' data(design)
#' data(covset)
#' Subject <- covset$ear
#' Time <- covset$time
#' Nboot <- 2  # for real data analysis, use Nboot at least 100
#' ncores <- 1 # for real data analysis and if the computer allows, increase ncores to save time
#' print.progress <- FALSE
#' saveboot <- FALSE
#' circadian <- TRUE
#' counts <- dat[1:3,]
#' C.matrix <- list()
#' # test for Line main effect
#' C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
#' # test for Time main effect
#' C.matrix[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
#' names(C.matrix) <- c("line2", "time")
#' TCout <- rmRNAseq:::TC_CAR1(counts, design, Subject, Time, C.matrix,
#' Nboot, ncores, print.progress, saveboot, circadian)
#' names(TCout)
#' TCout$pqvalue$pv
#' TCout$pqvalue$qv
#'
#' @export


TC_CAR1 <- function(counts, design, Subject, Time, C.matrix, Nboot = 100, ncores = 1, print.progress = FALSE, saveboot = FALSE, circadian = TRUE) {
  # message("Analyzing data using voomgls \n")

  v <- myvoom(counts = counts, design = design, lib.size = apply(counts, 2, stats::quantile,
                                                                 0.75), plot = FALSE)
  # Estimate new time point
  # message("Estimate new time point \n")
  if(circadian){
    tmOut <- TimeMin(v, Subject, Time, ncores)
  TimeMinOut <-tmOut$MinMaxTime
  ntp <- NewTimeEst(v, Subject, Time, TimeMinOut, ncores)
  NewTime <- Time
  TimeOther <- setdiff(unique(Time), TimeMinOut)
  NewTime[Time==TimeMinOut[1]] <- 1
  NewTime[Time==TimeMinOut[2]] <- 0
  NewTime[Time==TimeOther[1]] <- ntp[1]
  NewTime[Time==TimeOther[2]] <- ntp[2]
  Time <- NewTime}
  newlm0 <- voomgls_CAR1(v = v, Subject = Subject, Time = Time, ncores = ncores, C.matrix = C.matrix,
                         beta0 = NULL)


  newlm <- voomgls_Symm(v = v, Subject = Subject, Time = Time, ncores = ncores, C.matrix = C.matrix,
                        beta0 = NULL, print.progress)

  BetaMat <- data.matrix(newlm[grep("fixed.", names(newlm))])
  Sigma2Vec <- newlm$s2_shrunken
  RhoVec <- data.matrix(newlm[grep("rho", names(newlm))])
  WeightMat <- v$weights
  lib.size <- v$targets$lib.size

  # message("Analyzing bootstrap data using voomgls \n")
  bootres <- parallel::mclapply(1:Nboot, function(nrep){
    # if (print.progress)
    message("running bootstrap sample nrep = ", nrep, "\n")

    message("---------------------------------------------\n")

    simcounts <- sc_Symm(BetaMat = BetaMat, Sigma2Vec = Sigma2Vec, RhoVec = RhoVec, WeightMat = WeightMat,
                         lib.size = lib.size, design = design, Subject = Subject, Time = Time, nrep = nrep)
    # saveRDS(simcounts, file = paste0("simcounts_", nrep, ".rds"))
    v1 <- limma::voom(counts = simcounts, lib.size = apply(simcounts, 2, stats::quantile, 0.75),
                      design = design, plot = FALSE)
    bootnrep <- list(simcounts = simcounts, v1 = v1)
    bootlm <- voomgls_CAR1(v = v1, Subject = Subject, Time = Time, ncores = ncores, C.matrix = C.matrix,
                           beta0 = BetaMat, print.progress = print.progress)



    # if (print.progress)
    # message("finishing bootstrap sample nrep = ", nrep, "\n")

    cbind(simcounts, bootlm)
    # bootlm

  }, mc.cores = 1)
  # message("finishing analysis of all bootstrap samples \n")
  # node.error <- which(lapply(bootres, class) == "try-error")
  # if (length(node.error) != 0) {
  #   bootres <- res[lapply(bootres, class) != "try-error"]
  #   message("cores have error occured: ", node.error, "\n")
  # } else {
  #   # message("all cores work fine! \n")
  # }
  # calculating pvalue for all gene and all test in C.matrix components
  pvboot <- vapply(grep(pattern = "Ftest", names(newlm0), value = TRUE), function(x) {
    boot_test <- sapply(bootres, "[[", x)
    obs_test <- newlm0[, x]
    pv <- parallel::mclapply(1:length(obs_test), function(i) {
      temp <- boot_test >= obs_test[i]
      (sum(temp) + 1)/(length(temp) + 1)
    }, mc.cores = ncores)
    pv <- do.call("c", pv)
    pv
  }, FUN.VALUE = rep(1, nrow(counts)))
  # message("Finishing calculation of p-values and q values\n")
  pvboot <- data.frame(pvboot)
  names(pvboot) <- gsub("Ftest.", "", names(pvboot))
  qvboot <- data.frame(vapply(pvboot, function(x) jabes.q(x), FUN.VALUE = rep(1, nrow(counts))))
  pqvalue <- list(pv = pvboot, qv = qvboot)
  if(circadian){if (saveboot){res <- list(NewTime = NewTime, TimeMinOut = tmOut,
                            ori.res = list(v = v, newlmSymm = newlm,
                                           newlm = newlm0),
                            boot.res = bootres,
                            pqvalue = pqvalue)}else{
                              res <- list(NewTime = NewTime, TimeMinOut = tmOut,
                                          ori.res = list(v = v, newlmSymm = newlm,
                                                         newlm = newlm0),
                                          # boot.res = bootres,
                                          pqvalue = pqvalue)
                            }
  }else{
    if (saveboot){res <- list(ori.res = list(v = v, newlmSymm = newlm,
                                             newlm = newlm0),
                              boot.res = bootres,
                              pqvalue = pqvalue)}else{
                                res <- list(ori.res = list(v = v, newlmSymm = newlm,
                                                           newlm = newlm0),
                                            # boot.res = bootres,
                                            pqvalue = pqvalue)
                              }
  }
  # message("# p.05 of each test \n")
  # print(apply(pqvalue$pv <= 0.05, 2, sum))
  # message("\n # q.05 of each test \n")
  # print(apply(pqvalue$qv <= 0.05, 2, sum))
  return(res)
}
