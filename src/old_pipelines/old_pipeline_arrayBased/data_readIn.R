suppressMessages(lapply(c("dplyr","minfi"), library, character.only=T))

args <- commandArgs(trailingOnly=T)
job_idx <- as.numeric(args[1])
total_jobs <- as.numeric(args[2])

load("../data/targets.RData")

# Requires in the realm of 80Mb memory and 5 secs per sample
samples <- seq(from=floor((job_idx-1)/total_jobs*nrow(targets))+1,  # Sets "chunk" of sample indices to analyze
              to=floor(job_idx/total_jobs*nrow(targets)))

rgSet <- read.metharray.exp(targets=targets[samples,])

save("rgSet", file=paste0("../data/rgSet", job_idx, ".RData"))