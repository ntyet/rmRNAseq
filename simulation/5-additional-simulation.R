rm(list = ls())
pm1 <- proc.time()
library(dplyr)
library(splineTimeR)
library(ImpulseDE2)
library(rmRNAseq)

# This function conducts additional simulations using splineTimeR, ImpulseDE2,
# CAR1, and Symm



additional_sim <- function(scenario, effect, nrep, EE = 4000, DE = 1000, ncores = 20, name_dir_sim){
  covset <- rmRNAseq::covset

  # name_dir_sim <- paste0('./simulation/SimResults/Scenario_', scenario,
  #                        '_Effect_', effect)
  # name_dir_sim <- paste0('SimResults/Scenario_', scenario,
  #                        '_Effect_', effect)
  all_sim <- readRDS(paste0(name_dir_sim, '/All_sim_', nrep, '.rds'))
  #---results when applying CAR1----
  Ftest_CAR1 <- all_sim$resCAR1Symm$ori.res$newlm[, grep("Ftest.",
                                                         names(all_sim$resCAR1Symm$ori.res$newlm))]/(1*(effect==1) + 3*(effect ==2)) # 1 is the degree of freedom of the numerator
  pv_CAR1 <- 1-pf(q = Ftest_CAR1, df1 = 1*(effect==1) + 3*(effect ==2), df2 = 24)

  #---results when applying Symm----
  Ftest_Symm <- all_sim$resCAR1Symm$ori.res$newlmSymm[, grep("Ftest.",
                                                         names(all_sim$resCAR1Symm$ori.res$newlmSymm))]/(1*(effect==1) + 3*(effect ==2)) # 1 is the degree of freedom of the numerator
  pv_Symm <- 1-pf(q = Ftest_Symm, df1 = 1*(effect==1) + 3*(effect ==2), df2 = 24)

  if(effect == 1){
    #---results when applying splineTimeR---


    ##splineTimeR

    inData <- all_sim$simout$sim.counts
    colnames(inData) <- paste(covset$line, covset$ear, covset$timef,sep="_")

    design1 <- data.frame(row.names=colnames(inData),
                          "SampleName"=colnames(inData),
                          "Time"=covset$time,
                          "Treatment"=covset$line)
    phenoData <- new("AnnotatedDataFrame",data=design1)
    data <- ExpressionSet(assayData=as.matrix(inData),phenoData=phenoData)

    #---splineTimeR without log-transformed---

    diffExprs <-rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
                                              reference = "L",intercept = TRUE, voom_method = FALSE)

    #---splineTimeR with log-transformed---
    my_diffExprs <- rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
                                                  reference = "L",intercept = TRUE, voom_method = TRUE)

    #---results when applying ImpulseDE2---

    inData <- all_sim$simout$sim.counts
    colnames(inData) <- paste(covset$line, covset$ear, covset$timef,sep="_")
    row.names(inData) <- 1:nrow(inData)
    # specify experimental design
    design1 <- data.frame("Sample"=colnames(inData),
                          "Condition"=covset$line,
                          "Time"=covset$time,
                          row.names=colnames(inData))
    #--- for line effect---
    levels(design1$Condition) <- c("case", "control")
    # DEG analysis
    impulse_results <- runImpulseDE2(matCountData = inData,
                                     dfAnnotation =design1,
                                     boolCaseCtrl = TRUE,
                                     scaNProc = ncores,
                                     scaQThres = NULL,
                                     boolIdentifyTransients = TRUE,
                                     vecSizeFactorsExternal = NULL)
    pv <- data.frame(Symm = pv_Symm,
                     CAR1 = pv_CAR1,
                     splineTimeR = diffExprs$pv,
                     splineTimeRvoom = my_diffExprs$pv,
                     ImpulseDE2 = impulse_results$dfImpulseDE2Results$p)

    out_all <- list(splineTimeRout = diffExprs,
                    splineTimeRvoomout = my_diffExprs,
                    ImpulseDE2out = impulse_results)
    }else{
    #---results when applying ImpulseDE2---
    #--- for time effect---
    inData <- all_sim$simout$sim.counts
    colnames(inData) <- paste(covset$line, covset$ear, covset$timef,sep="_")
    row.names(inData) <- 1:nrow(inData)
    # specify experimental design
    design1 <- data.frame("Sample"=colnames(inData),
                          "Condition"=covset$line,
                          "Time"=covset$time,
                          row.names=colnames(inData))
    levels(design1$Condition) <- c("case", "case")
    # DEG analysis
    impulse_results <- runImpulseDE2(matCountData = inData,
                                     dfAnnotation =design1,
                                     boolCaseCtrl = FALSE,
                                     scaNProc = ncores,
                                     scaQThres = NULL,
                                     boolIdentifyTransients = TRUE,
                                     vecSizeFactorsExternal = NULL)
    pv <- data.frame(Symm = pv_Symm,
                     CAR1 = pv_CAR1,
                     ImpulseDE2 = impulse_results$dfImpulseDE2Results$p)
    out_all <- list(ImpulseDE2out = impulse_results)
  }
  oo <- vapply(pv, function(x)rmRNAseq:::pauc_out(x, EE = EE, DE = DE), FUN.VALUE = c(rep(1.0, 13)))
  oo <- reshape::melt(oo)
  out <- oo$value
  names(out) <- paste(oo$X1, oo$X2, sep  =".")
  Output_All_sim <- readRDS(paste0(name_dir_sim, '/Output_All_sim_', nrep, '.rds'))
  out <- c(Output_All_sim, out)
  saveRDS(out, file = paste0(name_dir_sim, "/Output_All_sim_", nrep, ".rds"))
  All_sim <- readRDS(paste0(name_dir_sim, '/All_sim_', nrep, '.rds'))
  out_all <- c(All_sim, out_all)
  saveRDS(out_all, file = paste0(name_dir_sim, "/All_sim_", nrep, ".rds"))
  out
  }



