#' A Wrap Function to analyze a Simulated Data - All Cases
#'
#' This function does the following: 1. Generating a simulated
#' counts data set consisting of EE  genes, and DE genes with
#' respect to some contrast from a ideal case, i.e., the counts
#' are generated from corCAR1 structure, or misspecified case, i.e., the counts
#' are generated from corSymm structure; analyzing this simulated data set
#' using  methods: \code{\link{TC_CAR1}},
#'  \code{\link{voomlimmaFit}}, \code{\link{edgeRFit}}
#' , \code{\link{DESeq2Fit}}.
#' @inheritParams voomgls_CAR1
#' @param EE number of EE genes
#' @param DE number of DE genes
#' @param nrep index of sim replicate
#' @param RFIanalysis the output from RFI RNA-seq dataset. In the ideal simulation case,
#' it is res from TC_CAR1, in the  misspecified case, it is res from TC_Symm
#' @param scenario either 2- 'Symm' or 1- 'CAR1'
#' @param saveboot TRUE or FALSE to save or not save bootstrap output
#' @param Subject a vector of subjects or experimental units.
#' @param Time a vector of time points.
#' @param C.matrix is a list of matrix Ci in testing H0:  Ci*beta = 0.
#' @param Nboot number of bootstrap replicates, default is 100.
#'@param print.progress \code{TRUE} or \code{FALSE}, printing the process or not.
#' @param name_dir_sim name of directory to contain the output and result of this function
#' @return R, V, FDR, PAUC, and AUC of all 7 methods (2 oracles with unshrunken and shrunken)
#' with FPR = 0.05, 0.10, 0.20 for each S, R, V, FDR, PAUC.
#' @importFrom reshape melt
#' @importFrom stats coef model.matrix quantile vcov
#' @examples
#' \donttest{
#' data(resSymm)
#' data(design)
#' data(covset)
#' C.matrix <- list()
#' # test for Line main effect
#' C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
#' names(C.matrix) <- c("line2")
#' EE <- 10; DE <- 5; ncores <- 1; Subject <- covset$ear;
#' Time <- covset$time; Nboot <- 2; nrep <- 1; name_dir_sim <- "Output/MisspecifiedCased";
#' print.progress=F
#' TC_Symm_scOut <- rmRNAseq:::TC_CAR1_sc(resSymm,  EE, DE,  C.matrix,
#' Subject, Time,  Nboot,  nrep, ncores, name_dir_sim , print.progress)
#' names(TC_Symm_scOut)
#' }
#'
TC_CAR1_sc <- function(RFIanalysis, scenario,  EE, DE,  C.matrix, Subject, Time,
                             Nboot,  nrep, ncores, name_dir_sim = NULL, print.progress=F,
                             saveboot = FALSE){
  if(!is.null(name_dir_sim)) dir.create(name_dir_sim, showWarnings = FALSE, recursive = TRUE)
  Cmatrixseed <- c(1, 2)
  names(Cmatrixseed) <- c("line2", "time")
  if( scenario == 2){
    res <- RFIanalysis$Symm
    set.seed(Cmatrixseed[names(C.matrix)]*10^7 + nrep)
    design <- res$ori.res$v$design
    pv <- res$pqvalue$pv[,names(C.matrix)]
    sortline <- sort(pv, index.return = T)
    DEline.ind <- sortline$ix[1:(length(pv)/2)]
    EEline.ind <- setdiff(1:length(pv), DEline.ind)
    coefbeta <- res$ori.res$newlm[,grep("fixed.", names(res$ori.res$newlm))]
    colnames(coefbeta) <- gsub("fixed.", "", colnames(coefbeta))
    coefbeta[EEline.ind, colnames(C.matrix[[1]])]  <- 0
    DE.sample <- sample(DEline.ind, size = DE, replace = F)
    EE.sample <- sample(EEline.ind, size = EE, replace = F)
    sim.sample <- c(EE.sample, DE.sample)
    names(sim.sample) <- rownames(res$ori.res$newlm)[sim.sample]
    BetaMat <- coefbeta[sim.sample,]
    Sigma2Vec <- res$ori.res$newlm$s2_shrunken[sim.sample]
    RhoVec <- data.matrix(res$ori.res$newlm[grep("rho", names(res$ori.res$newlm))])[sim.sample,]
    #WeightMat <- res$ori.res$v$weights[sim.sample,]
    lib.size <- res$ori.res$v$target$lib.size
    ## modified weights
    f <- res$ori.res$v$f
    fitted.values <- data.matrix(BetaMat) %*% t(design)
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    WeightMat <- w
    # plot(WeightMat, w, pch = .2)
    ## end modified weights
    simcounts <- sc_Symm(BetaMat, Sigma2Vec, RhoVec, WeightMat, lib.size, design, Subject, Time, nrep)
  }else{
    res <- RFIanalysis$CAR1
    set.seed(Cmatrixseed[names(C.matrix)]*10^6 + nrep)
    design <- res$ori.res$v$design
    pv <- res$pqvalue$pv[,names(C.matrix)]
    qv <- res$pqvalue$qv[,names(C.matrix)]
    # sum(qline <= .05)
    sortline <- sort(pv, index.return = T)
    DEline.ind <- sortline$ix[1:(length(pv)/2)]
    EEline.ind <- setdiff(1:length(pv), DEline.ind)
    coefbeta <- res$ori.res$newlm[,grep("fixed.", names(res$ori.res$newlm))]
    colnames(coefbeta) <- gsub("fixed.", "", colnames(coefbeta))
    coefbeta[EEline.ind, colnames(C.matrix[[1]])]  <- 0
    DE.sample <- sample(DEline.ind, size = DE, replace = F)
    EE.sample <- sample(EEline.ind, size = EE, replace = F)
    sim.sample <- c(EE.sample, DE.sample)
    names(sim.sample) <- rownames(res$ori.res$v$E)[sim.sample]
    BetaMat <- coefbeta[sim.sample,]
    Sigma2Vec <- res$ori.res$newlm$s2_shrunken[sim.sample]
    RhoVec <- res$ori.res$newlm$rho[sim.sample]
    #WeightMat <- res$ori.res$v$weights[sim.sample,]
    lib.size <- res$ori.res$v$target$lib.size
    ## modified weights
    f <- res$ori.res$v$f
    fitted.values <- data.matrix(BetaMat) %*% t(design)
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    WeightMat <- w
    ## end modified weights

    simcounts <- sc_CAR1(BetaMat, Sigma2Vec, RhoVec, WeightMat, lib.size, design, Subject, Time, nrep)

  }
  simout <-  list(sim.counts = simcounts, sim.sample = sim.sample)

  ## Reestimate time points
  cat("--------------------------------------\n")

  cat("Analyzing simcounts using TC-CAR1 \n")

  resCAR1Symm <- TC_CAR1(counts = simcounts, design = design,
                              Subject = Subject, Time = Time,  C.matrix,
                              Nboot = Nboot,  ncores = ncores, print.progress, saveboot)
  cat("--------------------------------------\n")

  cat("Analyzing simcounts using voomlimma \n")

  voomlimmaout <- voomlimmaFit(simcounts, design, names(C.matrix))
  cat("--------------------------------------\n")

  cat("Analyzing simcounts using edgeR \n")

  edgeRout <- edgeRFit(simcounts, design, names(C.matrix))
  cat("--------------------------------------\n")

  cat("Analyzing simcounts using DESeq2 \n")

  DESeq2out <- suppressMessages(DESeq2Fit(simcounts, design, names(C.matrix)))
  cat("--------------------------------------\n")

  ressim <- list(simout = simout,
    resCAR1Symm=resCAR1Symm,
    voomlimmaout = voomlimmaout,
    edgeRout = edgeRout,
    DESeq2out = DESeq2out
  )
  if(!is.null(name_dir_sim))saveRDS(ressim, file = paste0(name_dir_sim, "/All_sim_", nrep, ".rds"))

  pv <- data.frame(bootCAR1Symm = ressim$resCAR1Symm$pqvalue$pv[,names(C.matrix)],
                   voomlimma = ressim$voomlimmaout$pv[, names(C.matrix)],
                   edgeR = ressim$edgeRout$pv[, names(C.matrix)],
                   DESeq2 = ressim$DESeq2out$pv[, names(C.matrix)])

  oo <- vapply(pv, function(x)pauc_out(x, EE = EE, DE = DE), FUN.VALUE = c(rep(1.0, 13)))
  oo <- melt(oo)
  out <- oo$value
  names(out) <- paste(oo$X1, oo$X2, sep  =".")
  if(!is.null(name_dir_sim))saveRDS(out, file = paste0(name_dir_sim, "/Output_All_sim_", nrep, ".rds"))
  out
}

