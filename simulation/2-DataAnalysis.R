rm(list = ls())
pm1 <- proc.time()
library(splineTimeR)
library(ImpulseDE2)
library(rmRNAseq)
data(design)
alldata <- readRDS("simulation/data/rfidata.rds")
counts <- alldata$counts
Subject <- alldata$covset$ear
Time <- alldata$covset$time
Nboot <- 100
ncores <- 20
print.progress <- F
C.matrix0 <- list()
C.matrix0[[1]] <- limma::makeContrasts(line2, levels = design)
C.matrix0[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
names(C.matrix0) <- c('line2', 'time')

pm <- proc.time()
bootCAR1 <- rmRNAseq::TC_CAR1(counts = counts, design = design, Subject = Subject,
                                Time = Time, C.matrix = C.matrix0, Nboot = Nboot,
                                ncores = ncores, print.progress = print.progress)
proc.time() - pm

pm <- proc.time()
bootSymm <- rmRNAseq::TC_Symm(counts = counts, design = design, Subject = Subject,
                                Time = Time, C.matrix = C.matrix0, Nboot = Nboot,
                                ncores = ncores,print.progress = print.progress)
proc.time() - pm

pm <- proc.time()
voomout <- plyr::llply(names(C.matrix0),
            function(effect)rmRNAseq:::voomlimmaFit(counts = counts,
                                                      design = design,
                                                      Effect = effect))
proc.time() - pm
pm <- proc.time()
edgeRout <- plyr::llply(names(C.matrix0),
                       function(effect)rmRNAseq:::edgeRFit(counts = counts,
                                                                 design = design,
                                                                 Effect = effect))
proc.time() - pm
pm <- proc.time()
DESeq2out <- plyr::llply(names(C.matrix0),
                        function(effect)rmRNAseq:::DESeq2Fit(counts = counts,
                                                              design = design,
                                                              Effect = effect))
proc.time() - pm
names(voomout) <- names(edgeRout) <- names(DESeq2out) <- names(C.matrix0)



##splineTimeR

inData <- alldata$counts
colnames(inData) <- paste(covset$line, covset$ear, covset$timef,sep="_")

design1 <- data.frame(row.names=colnames(inData),
                      "SampleName"=colnames(inData),
                      "Time"=covset$time,
                      "Treatment"=covset$line)
phenoData <- new("AnnotatedDataFrame",data=design1)
data <- ExpressionSet(assayData=as.matrix(inData),phenoData=phenoData)

#---splineTimeR without log-transformed---
pm <- proc.time()
diffExprs <-rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
                                          reference = "L",intercept = TRUE, voom_method = FALSE)
proc.time() - pm

#---splineTimeR with log-transformed---
pm <- proc.time()
my_diffExprs <- rmRNAseq:::my_splineDiffExprs(eSetObject = data, df = 3,cutoff.adj.pVal = 1,
                                              reference = "L",intercept = TRUE, voom_method = TRUE)

proc.time() - pm

#---ImpulseDE2---

pm <- proc.time()
inData <- as.matrix(alldata$counts)
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
impulse_results_line <- runImpulseDE2(matCountData = inData,
                                 dfAnnotation =design1,
                                 boolCaseCtrl = TRUE,
                                 scaNProc = 4,
                                 scaQThres = NULL,
                                 boolIdentifyTransients = TRUE,
                                 vecSizeFactorsExternal = NULL)
proc.time() - pm
#--- for time effect---
inData <- as.matrix(alldata$counts)
colnames(inData) <- paste(covset$line, covset$ear, covset$timef,sep="_")
row.names(inData) <- 1:nrow(inData)
# specify experimental design
design1 <- data.frame("Sample"=colnames(inData),
                      "Condition"=covset$line,
                      "Time"=covset$time,
                      row.names=colnames(inData))
levels(design1$Condition) <- c("case", "case")
# DEG analysis
impulse_results_time <- runImpulseDE2(matCountData = inData,
                                 dfAnnotation =design1,
                                 boolCaseCtrl = FALSE,
                                 scaNProc = 4,
                                 scaQThres = NULL,
                                 boolIdentifyTransients = TRUE,
                                 vecSizeFactorsExternal = NULL)

res <- list(CAR1 = bootCAR1,
            Symm = bootSymm,
            voom = voomout,
            edgeR = edgeRout,
            DESeq2 = DESeq2out,
            splineTimeR = diffExprs,
            splineTimeRvoom = my_diffExprs,
            ImpulseDE2 = list(line = impulse_results_line, time = impulse_results_time))

# res <- list(splineTimeR = diffExprs,
#             splineTimeRvoom = my_diffExprs,
#             ImpulseDE2 = list(line = impulse_results_line, time = impulse_results_time))

# RFIanalysis <- readRDS("simulation/DataAnalysis/RFIanalysis.rds")
# res <- c(RFIanalysis, res)
saveRDS(res, file = paste0("simulation/DataAnalysis/RFIanalysis.rds"))
proc.time() - pm1

