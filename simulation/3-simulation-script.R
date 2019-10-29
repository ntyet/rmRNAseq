# setwd("./simulation/")
EE <- 4000; DE <- 1000; Nboot <- 100; ncores <- 20; ns <- 10; Seqsim <-0:4;  walltime <- "2-12:00:00"

Effects <- c("line2", "time")
listEffects <- list("line2", "time2, time6, time24")
names(listEffects) <- Effects
Scenario <- c("CAR1", "Symm")
for (seqsim in Seqsim){
  for (scenario in 2:1){
    for (effect in 2:1){
    filetext <- paste0("rm(list = ls())
pm1 <- proc.time()
library(rmRNAseq)
C.matrix0 <- list()
C.matrix0[[1]] <- limma::makeContrasts(line2, levels = design)
C.matrix0[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
names(C.matrix0) <- c('line2', 'time')
C.matrix <- C.matrix0[", effect, "]
EE <- ", EE, "; DE <- ", DE, "; Subject <- covset$ear; Time <- covset$time; Nboot <- ", Nboot, "; ncores <- ", ncores, ";  nSim <- ", seqsim*ns, "+1:", ns, "; print.progress = F; saveboot = T
RFIanalysis <- readRDS('DataAnalysis/RFIanalysis.rds')
res <- RFIanalysis$",Scenario[scenario], "
cat('scenario = ', ", scenario,", 'effect = ', ", effect, ", 'nsim =',  min(nSim), '--',  max(nSim), '\\n')
name_dir_sim <- 'SimResults/Scenario_", scenario, "_Effect_", effect, "'
# ##---------------------------------------------------------------------
output <- lapply(nSim, function(nrep){
cat('nrep = ', nrep, '\\n')
rmRNAseq:::TC_CAR1_sc(RFIanalysis, ", scenario, ",EE, DE, C.matrix, Subject, Time, Nboot, nrep, ncores,
name_dir_sim, print.progress, saveboot)
})
rm(output)
# output <- do.call('rbind', output)
# saveRDS(output, file = paste0(name_dir_sim, '/All_output_', min(nSim), '_', max(nSim), '.rds'))
proc.time() - pm1
")

filebash <- paste0("#!/bin/bash -l
#Submit this script with: sbatch thefilename
#SBATCH -t ", walltime, "   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n ", ncores, "   # total number of processor cores in this job
#SBATCH -J 'CAR1Symm-sno-",scenario, "-eft-",effect,  "-nsm-", seqsim*ns+1, "-", seqsim*ns+ns, "'   # job name
#SBATCH --mail-user=tienyettoan@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Load R module
enable_lmod
module load icc
module load libstdcxx/4
module load R
R CMD BATCH CAR1Symm-sno-",scenario, "-eft-",effect,  "-nsm-", seqsim*ns+1, "-", seqsim*ns+ns, ".R
exit

")


cat(filetext,file=paste0('CAR1Symm-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.R'))
cat(filebash,file=paste0('CAR1Symm-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
system(paste0('sbatch CAR1Symm-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
print(paste0('sbatch CAR1Symm-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
    }
  }
}

