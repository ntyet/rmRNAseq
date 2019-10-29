# This function conducts additional simulations using splineTimeR, ImpulseDE2,
# CAR1, and Symm

EE <- 4000; DE <- 1000; ncores <- 4;
ns <- 5; Seqsim <-0:9;  walltime <- "3-12:00:00"
for(seqsim in Seqsim){
  for(scenario in 1:2){
    for(effect in 1:2){
    filetext <- paste0("rm(list = ls())
library(rmRNAseq)
library(splineTimeR)
library(ImpulseDE2)
source('5-additional-simulation.R')
EE <- ", EE, "; DE <- ", DE, "; ncores <- ", ncores, "; nSim <- ", seqsim*ns, "+1:", ns,";
cat('scenario = ', ", scenario,", 'effect = ', ", effect, ", 'nsim =',  min(nSim), '--',  max(nSim), '\\n')

name_dir_sim <- 'SimResults/Scenario_", scenario, "_Effect_", effect, "'

# ##---------------------------------------------------------------------

output <- lapply(nSim, function(nrep){
cat('nrep = ', nrep, '\\n')
additional_sim(", scenario, ",", effect, ",","nrep, EE, DE, ncores, name_dir_sim)
})
proc.time() - pm1
")

filebash <- paste0("#!/bin/bash -l
#Submit this script with: sbatch thefilename
                       #SBATCH -t ", walltime, "   # walltime
                       #SBATCH -N 1   # number of nodes in this job
                       #SBATCH -n ", ncores, "   # total number of processor cores in this job
                       #SBATCH -J 'add-sim-sno-",scenario, "-eft-",effect,  "-nsm-", seqsim*ns+1, "-", seqsim*ns+ns, "'   # job name
                       #SBATCH --mail-user=tienyettoan@gmail.com   # email address
                       #SBATCH --mail-type=BEGIN
                       #SBATCH --mail-type=END
                       #SBATCH --mail-type=FAIL

                       # Load R module
                       enable_lmod
                       module load icc
                       module load libstdcxx/4
                       module load R
                       R CMD BATCH add-sim-sno-",scenario, "-eft-",effect,  "-nsm-", seqsim*ns+1, "-", seqsim*ns+ns, ".R
                       exit

                       ")

cat(filetext,file=paste0('add-sim-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.R'))
cat(filebash,file=paste0('add-sim-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
system(paste0('sbatch add-sim-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
print(paste0('sbatch add-sim-sno-', scenario, '-eft-',effect, '-nsm-', seqsim*ns+1, '-', seqsim*ns+ns,'.sh'))
    }
  }
  }
