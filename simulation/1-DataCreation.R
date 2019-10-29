# BEGIN DATA PREPARATION --------------------------------------------------------
library(openxlsx)
library(plyr)
library(dplyr)
#read meta data containing covariates
metadata <- read.xlsx("simulation/data-raw/ISS.RNA-seq.final.metadata.for.Yet v2.xlsx")
covset <- metadata%>%
  dplyr::arrange(Ear.Tag, Time)%>% # reorder by incresing time, and ear tag
  dplyr::filter(Project_II_data ==2)%>%
  dplyr::mutate(sampleid = as.character(Sample.ID),
                batch = as.factor(as.character(Sequencin_batch)),
                rnadate = as.factor(as.character(RNA_prep_date)),
                ear = as.numeric(as.character(Ear.Tag)),
                BW2 = ifelse(ear==8, 65.6,
                             ifelse(ear ==16, 65.8,
                                    ifelse(ear ==19, 70.8,
                                           ifelse(ear==41, 69.6,
                                                  ifelse(ear ==50, 65.6,
                                                         ifelse(ear==139, 68.8,
                                                                ifelse(ear ==573, 61.0, 62.4))))))),
                time = as.numeric(as.character(Time)),
                timef = as.factor(Time),
                line = as.factor(as.character(RFI_line)),
                rina = as.numeric(as.character(post_RIN)),
                neut = as.numeric(as.character(Neutrophil)),
                lymp = as.numeric(as.character(Lymphocyte)),
                mono = as.numeric(as.character(Monocyte)),
                baso = as.numeric(as.character(Basophils)),
                eosi = as.numeric(as.character(Eosinophil)),
                plat = as.numeric(as.character(Platelet_Auto)),
                rbc = as.numeric(as.character(RBC)),
                wbc = as.numeric(as.character(WBC)),
                total_cell = plat + wbc + rbc,
                lrina = log(rina),
                lBW2 = log(BW2),
                tot_cell = wbc + rbc + plat,
                pneut = neut/tot_cell,
                plymp = lymp/tot_cell,
                pbaso = baso/tot_cell,
                peosi = eosi/tot_cell,
                pmono = mono/tot_cell,
                pplat = plat/tot_cell,
                prbc = rbc/tot_cell,
                lpneut = log(neut/tot_cell),
                lplymp = log(lymp/tot_cell),
                lpmono = log(mono/tot_cell),
                lpbaso = log(baso/tot_cell),
                lpeosi = log(eosi/tot_cell),
                lprbc = log(rbc/tot_cell),
                lpplat = log(plat/tot_cell),
                pneutW = neut/wbc,
                plympW = lymp/wbc,
                pbasoW = baso/wbc,
                peosiW = eosi/wbc,
                pmonoW = mono/wbc,
                lpneutW = log(neut/wbc),
                lplympW = log(lymp/wbc),
                lpmonoW = log(mono/wbc),
                lpbasoW = log(baso/wbc),
                lpeosiW = log(eosi/wbc)
  )%>%
  dplyr::select(time, timef, line, ear, sampleid, batch, rnadate, BW2, lBW2, rina, neut, lymp, mono, baso, eosi, rbc, plat, lpneut, lrina,
                lplymp, lpmono, lpbaso, lpeosi, lprbc, lpplat, tot_cell, plymp, pneut,pmono, pbaso, peosi, pplat, prbc,
                plympW, pneutW,pmonoW, pbasoW, peosiW, lplympW, lpneutW,lpmonoW, lpbasoW, lpeosiW, wbc)


covset <- within(covset, {
  lpbaso[which(!is.finite(lpbaso))] <- log((baso[which(!is.finite(lpbaso))]+.01/2)/tot_cell[which(!is.finite(lpbaso))])
  lpbasoW[which(!is.finite(lpbasoW))] <- log((baso[which(!is.finite(lpbasoW))]+.01/2)/wbc[which(!is.finite(lpbasoW))])
})


# design mat
CC <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0,
               0, 1, 0, 0, 0, 1/4, 1/4, 1/4, # line effect
               0, 0, 1, 0, 0, 1/2, 0, 0, # time effect T2 - T0
               0, 0, 0, 1, 0, 0, 1/2, 0, # time effect T6 - T0
               0, 0, 0, 0, 1, 0, 0, 1/2, # time effect T24 - T0
               0, 0, 0, 0, 0, 1, 0, 0, # line*time effect (T2L- T0L) - (T2H - T0H)
               0, 0, 0, 0, 0, 0, 1, 0, # line*time effect (T6L- T0L) - (T6H - T0H)
               0, 0, 0, 0, 0, 0, 0, 1), ncol = 8, byrow = T)  # line*time effect (T24L- T0L) - (T24H - T0H)
X <- model.matrix(~line*timef, data = covset)
Y <- X%*%solve(CC)
colnames(Y) <- c("Intercept", "line2", "time2", "time6", "time24",
                 "linetime2", "linetime6", "linetime24")
covset <- cbind(covset, Y)
covset$times <-covset$time/max(covset$time)
design0 <- model.matrix(~line2+time2 + time6 + time24 +
                          linetime2 + linetime6 + linetime24, data =covset)


colnames(design0)[1] <- "Intercept"
# rna-seq read count data set
dat <- read.table("simulation/data-raw/ISS.RNA-seq.count.final.table.txt",
                  header = T)

dat<- dat%>%
  arrange(Geneid)
rownames(dat) <- dat[, "Geneid"]
dat <- dat[, covset$sampleid]
dat<-dat[rowMeans(dat) > 8 &rowSums(dat!=0) > 3, ]
rownames(dat) <- NULL
colnames(dat) <- NULL
covset <- covset[c("time", "timef", "line", "ear", "line2", "time2", "time6", "time24", "linetime2", "linetime6", "linetime24")]
design <- design0

rfidata <- list(counts = dat, design0= design0, covset = covset)
saveRDS(rfidata, file = "simulation/data/rfidata.rds")

