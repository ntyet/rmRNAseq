###Sample Correlation plots-------
library(plyr)
library(tidyverse)
library(limma)
library(dplyr)
RFIanalysis <- readRDS("simulation/DataAnalysis/RFIanalysis.rds")
rho <- data.frame(rho1 = RFIanalysis$Symm$ori.res$newlm$rho1,
                  rho2 = RFIanalysis$Symm$ori.res$newlm$rho2,
                  rho3 = RFIanalysis$Symm$ori.res$newlm$rho3,
                  rho4 = RFIanalysis$Symm$ori.res$newlm$rho4,
                  rho5 = RFIanalysis$Symm$ori.res$newlm$rho5,
                  rho6 = RFIanalysis$Symm$ori.res$newlm$rho6) # (0, 2), (0, 6), (0, 24), (2, 6), (2, 24), (6, 26)

rho2 <- rho%>%gather(., key = "pairs", value = "correlation")%>%
  mutate(pairs = replace(pairs, pairs =="rho1", '(0,2)'))%>%
  mutate(pairs = replace(pairs, pairs =="rho2", '(0,6)'))%>%
  mutate(pairs = replace(pairs, pairs =="rho3", '(0,24)'))%>%
  mutate(pairs = replace(pairs, pairs =="rho4", '(2,6)'))%>%
  mutate(pairs = replace(pairs, pairs =="rho5", '(2,24)'))%>%
  mutate(pairs = replace(pairs, pairs =="rho6", '(6,24)'))

ggplot(data = rho2, mapping = aes( x = reorder(pairs, correlation, FUN = median), y = correlation))+
  geom_boxplot()+
  xlab("Time Pair")+
  ylab("Estimated Correlation")

ggsave("bioinformatics-rnaseq-repeated/SampleCorreltion_AllGenes.pdf", width = 6, height = 3)
ggsave("bioinformatics-rnaseq-repeated/SampleCorreltion_AllGenes.eps", width = 6, height = 3)
ggsave("simulation/figure/SampleCorreltion_AllGenes.pdf", width = 6, height = 3)
ggsave("simulation/figure/SampleCorreltion_AllGenes.eps", width = 6, height = 3)


{ sink("/dev/null"); str(RFIanalysis$ImpulseDE2$time); sink(); }


ImpulseDE2_time_pv <- RFIanalysis$ImpulseDE2$time$dfImpulseDE2Results$p
ImpulseDE2_line_pv <- RFIanalysis$ImpulseDE2$line$dfImpulseDE2Results$p
ImpulseDE2_time_qv <- rmRNAseq:::jabes.q(ImpulseDE2_time_pv)
ImpulseDE2_line_qv <- rmRNAseq:::jabes.q(ImpulseDE2_line_pv)

Line <- cbind(rmRNAseq = RFIanalysis$CAR1$pqvalue$qv$line2 <= .05,
              voomlimma = RFIanalysis$voom$line2$qv[,1] <= .05,
              edgeR = RFIanalysis$edgeR$line2$qv[,1] <= .05,
              DESeq2 = RFIanalysis$DESeq2$line2$qv[,1] <= .05,
              splineTimeR = RFIanalysis$splineTimeR$qv <=.05,
              ImpulseDE2 = ImpulseDE2_line_qv <= .05)

Time <- cbind(rmRNAseq = RFIanalysis$CAR1$pqvalue$qv$time <= .05,
              voomlimma = RFIanalysis$voom$time$qv[,1] <= .05,
              edgeR = RFIanalysis$edgeR$time$qv[,1] <= .05,
              DESeq2 = RFIanalysis$DESeq2$time$qv[,1] <= .05,
              ImpulseDE2 = ImpulseDE2_time_qv <= .05)

###Bar Plot number of DEGs-------

DEGs <- data.frame(DEGs = c(colSums(Line), colSums(Time)),
                   Method = factor(c(colnames(Line), colnames(Time)),
                                   levels = c("rmRNAseq", "voomlimma",
                                              "edgeR", "DESeq2",
                                              "splineTimeR", "ImpulseDE2")),
                   Effect = c(rep("Line", ncol(Line)), rep("Time", ncol(Time))))


ggplot(data = DEGs, mapping = aes(x = Method, y = DEGs)) +
  geom_bar(stat = "identity")+
  facet_wrap(~Effect, scale = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Number of DE genes")
theme(legend.position="none") #+
#scale_y_log10()

ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.pdf", width = 6, height = 3)
ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.png", width = 6, height = 3)
ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.eps", width = 6, height = 3)

ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.pdf", width = 6, height = 3)
ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.png", width = 6, height = 3)
ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.eps", width = 6, height = 3)

ggsave("simulation/figure/BarDiagramRFI_Time024.pdf", width = 6, height = 3)
ggsave("simulation/figure/BarDiagramRFI_Time024.png", width = 6, height = 3)
ggsave("simulation/figure/BarDiagramRFI_Time024.eps", width = 6, height = 3)


ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.pdf", width = 10, height = 7, units = "cm")
ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.png", width = 10, height = 7, units = "cm")
ggsave("bioinformatics-submitted-rv/BarDiagramRFI_Time024.eps", width = 10, height = 7, units = "cm")

ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.pdf", width = 10, height = 7, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.png", width = 10, height = 7, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/BarDiagramRFI_Time024.eps", width = 10, height = 7, units = "cm")

ggsave("simulation/figure/BarDiagramRFI_Time024.pdf", width = 10, height = 7, units = "cm")
ggsave("simulation/figure/BarDiagramRFI_Time024.png", width = 10, height = 7, units = "cm")
ggsave("simulation/figure/BarDiagramRFI_Time024.eps", width = 10, height = 7, units = "cm")

###Venn Diagrams-------

# pdf("bioinformatics-rnaseq-repeated/VennDiagramRFI_Time024.pdf", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()
#
#
# postscript("bioinformatics-rnaseq-repeated/VennDiagramRFI_Time024.eps", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()
#
# pdf("bioinformatics-rnaseq-repeated/figure/VennDiagramRFI_Time024.pdf", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()
#
#
# postscript("bioinformatics-rnaseq-repeated/figure/VennDiagramRFI_Time024.eps", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()
#
# pdf("simulation/figure/VennDiagramRFI_Time024.pdf", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()
#
#
# postscript("simulation/figure/VennDiagramRFI_Time024.eps", width = 11, height = 4)
# par(mfrow = c(1, 2))
# vennDiagram(vennCounts(Line))
# title("Line", adj = .5, line = 2, cex.main=2)
# vennDiagram(vennCounts(Time))
# title("Time", adj = .5, line = 2, cex.main=2)
# dev.off()



###Simulation Results-------

dirpath <- list.dirs(path = "simulation/SimResults")
PAUC <- NTP <- FDR <- NULL
for(dp in dirpath[-1]){ # dp <- dirpath[5]
  pout <- list.files(path = dp, pattern = "Output_All_sim", include.dirs = TRUE)
  pout <- paste0("Output_All_sim_", 1:50, ".rds")
  out <- plyr::laply(pout, function(pouti){ # pouti <- pout[5]
    fname <- paste0(dp,'/', pouti)
    fout <- readRDS(fname)
    fout
  })
  out1 <- data.frame(out,
                     Scenario = strsplit(dp, split = c("_"))[[1]][2],
                     Effect = strsplit(dp, split = c("_"))[[1]][4])
  colnames(out1) <- gsub("^auc.", "pauc1.", colnames(out1))
  fpr <- c(".05", ".10", ".20")
  # NTP0 <- llply(fpr, function(x){
  #   out3 <- out1%>%
  #     gather(key = "Method", value = "S",
  #            starts_with(paste0("S", x, collapse = "")) )%>%
  #     mutate(FPR = x,
  #            Method = gsub(pattern = paste0("S", x, ".", collapse = ""),
  #                          "", Method))%>%
  #     select(Method, FPR, Scenario, Effect, S)
  # })
  # NTP0 <- do.call("rbind", NTP0)


  R0 <- llply(fpr, function(x){ # x <- fpr[1]
    out3 <- out1%>%
      gather(key = "Method", value = "R",
             starts_with(paste0("R", x, collapse = "")) )%>%
      mutate(FPR = x,
             Method = gsub(pattern = paste0("R.", x, ".", collapse = ""),
                           "", Method))%>%
      dplyr::select(Method, FPR, Scenario, Effect, R)
  })
  R0 <- do.call("rbind", R0)

  V0 <- llply(fpr, function(x){
    out3 <- out1%>%
      gather(key = "Method", value = "V",
             starts_with(paste0("V", x, collapse = "")) )%>%
      mutate(FPR = x,
             Method = gsub(pattern = paste0("V.", x, ".", collapse = ""),
                           "", Method))%>%
      dplyr::select(Method, FPR, Scenario, Effect, V)
  })
  V0 <- do.call("rbind", V0)
  NTP0 <- R0
  NTP0$R <- R0$R - V0$V
  colnames(NTP0)[5] <- "S"

  FDR0 <- llply(fpr, function(x){
    out3 <- out1%>%
      gather(key = "Method", value = "FDR",
             starts_with(paste0("FDR", x, collapse = "")) )%>%
      mutate(FPR = x,
             Method = gsub(pattern = paste0("FDR", x, ".", collapse = ""),
                           "", Method))%>%
      dplyr::select(Method, FPR, Scenario, Effect, FDR)
  })
  FDR0 <- do.call("rbind", FDR0)
  fpr <- c(".05", ".10", ".20")

  PAUC0 <- llply(fpr, function(x){
    out3 <- out1%>%
      gather(key = "Method", value = "PAUC",
             starts_with(paste0("pauc", x, collapse = "")) )%>%
      mutate(FPR = x,
             Method = gsub(pattern = paste0("pauc", x, ".", collapse = ""),
                           "", Method))%>%
      dplyr::select(Method, FPR, Scenario, Effect, PAUC)
  })
  PAUC0 <- do.call("rbind", PAUC0)
  PAUC <- rbind(PAUC, PAUC0)
  NTP <- rbind(NTP, NTP0)
  FDR <- rbind(FDR, FDR0)
}


alloutput <- data.frame(FDR, NTP = NTP$S, PAUC = PAUC$PAUC)
levels(alloutput$Scenario) <- c("corCAR1", "corSymm")
levels(alloutput$Effect) <- c("Line", "Time", "Interaction")
d1 <- alloutput%>%mutate(Method = ifelse(Method == "bootCAR1Symm",  "rmRNAseq", Method))
d1$Method <- factor(d1$Method, levels = c("rmRNAseq",  "CAR1","Symm", "voomlimma", "edgeR", "DESeq2",
                                          "splineTimeRvoom",  "splineTimeR", "ImpulseDE2"))

# d1%>%dplyr::filter(FPR == ".05" & Effect == "Time" & Scenario == "corSymm")%>% dplyr::group_by(Method, FPR, Scenario, Effect)%>%dplyr::summarise(Ave = mean(NTP))

p1 <- ggplot(data = subset(d1, #Method %in% c("rmRNAseq", "voom", "edgeR", "DESeq2")&
                           (FPR == ".05")&
                             (Effect %in%  c("Line", "Time"))))+
  geom_boxplot(mapping = aes(y = FDR, x= Method))+
  facet_grid(Effect~Scenario, scale = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none")+
  geom_hline(aes(yintercept=.05))

ggsave("bioinformatics-submitted-rv/FDR_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/FDR_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/FDR_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("bioinformatics-rnaseq-repeated/FDR_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/FDR_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/FDR_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("simulation/figure/FDR_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/FDR_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/FDR_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")


p2 <- ggplot(data = subset(d1, #Method %in% c("rmRNAseq", "voom", "edgeR", "DESeq2")&
                           (FPR == ".05")&
                             (Effect  %in%  c("Line", "Time"))))+
  geom_boxplot(mapping = aes(y = PAUC, x= Method))+
  facet_grid(Effect~Scenario, scale = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none")

ggsave("bioinformatics-submitted-rv/PAUC_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/PAUC_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/PAUC_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("bioinformatics-rnaseq-repeated/PAUC_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/PAUC_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/PAUC_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("simulation/figure/PAUC_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/PAUC_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/PAUC_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")


p3 <- ggplot(data = subset(d1, #Method %in% c("rmRNAseq", "voom", "edgeR", "DESeq2")&
                           (FPR == ".05")&
                             (Effect  %in%  c("Line", "Time"))))+
  geom_boxplot(mapping = aes(y = NTP, x= Method))+
  facet_grid(Effect~Scenario, scale = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none")

ggsave("bioinformatics-submitted-rv/NTP_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/NTP_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-submitted-rv/NTP_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("bioinformatics-rnaseq-repeated/NTP_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/NTP_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("bioinformatics-rnaseq-repeated/NTP_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")

ggsave("simulation/figure/NTP_SimPlot_Time024.pdf", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/NTP_SimPlot_Time024.png", width = 10, height = 10, units = "cm")
ggsave("simulation/figure/NTP_SimPlot_Time024.eps", width = 10, height = 10, units = "cm")
