suppressMessages(silent <- lapply(c("tidyverse","minfi","wateRmelon","RefFreeEWAS","FlowSorted.Blood.450k"),
                        library, character.only=T))

## Sample sheet
targets <- read.metharray.sheet(base="../data/", pattern="^sample_sheet4.csv$") %>%
  dplyr::filter(Basename!="character(0)",  # Seems to appear in read.metharray.sheet output when a channel is missing
                Sample_Name!=21084)  # This sample is missing red channel .idat

## Get index and total number of jobs in array as well as # cores available from script arguments
args <- commandArgs(trailingOnly=T)
job_idx <- as.numeric(args[1])
total_jobs <- as.numeric(args[2])
num_cores <- as.numeric(args[3])

## Read in .idat files
# Requires in the realm of 75Mb memory and 5 secs per sample
samples <- seq(from=floor((job_idx-1)/total_jobs*nrow(targets))+1,  # Sets "chunk" of sample indices to analyze
               to=floor(job_idx/total_jobs*nrow(targets)))
print(paste("Samples to read:", min(samples), "to", max(samples)))

rgSet <- read.metharray.exp(targets=targets[samples,])
save("rgSet", file=paste0("../int/rgSet_", job_idx, ".RData"))  # Save as .RData with associated index no.
print("Samples read and saved.")
p3 <- proc.time()

## Generate probe detection p-values
detP <- detectionP(rgSet)
print("Detection-based QC run. Mem used so far:")
pryr::mem_used()
p4 <- proc.time()
print(p4-p3)

## Initial preprocessing with Noob (background adj. & dye bias)
mSet <- preprocessNoob(rgSet)
stopifnot(dim(detP)==dim(mSet))
"Preprocessing done."
rm(rgSet)

## Filter out bad samples and unwanted probes
keep_samples <- colSums(detP>0.01)/nrow(detP) < 0.1  # Don't want >10% of probes whose detection failed at p=0.01
detP <- detP[,keep_samples]  # Don't look at probes until bad samples have been removed
bad_probes <- rownames(detP)[rowSums(detP>0.01)/ncol(detP) > 0.01]  # Save these and filter out later
save("bad_probes", file=paste0("../int/bad_probes_", job_idx, ".RData"))
mSet.filt <- mSet[,keep_samples]  # Remove low-quality samples
print(paste("Removed", sum(!keep_samples), "samples"))
mSet.filt <- dropMethylationLoci(mSet.filt)  # Drop SNP probes and CpH probes
rm(detP, mSet)
print("QC complete. Mem used so far:")
pryr::mem_used()
p5 <- proc.time()
print(p5-p4)

## Run BMIQ within-sample normalization on each sample
# Requires approximately 0.5 Gb memory and 4 seconds per sample
betas.filt.norm_list <- mclapply(1:ncol(mSet.filt), function(x) BMIQ(mSet.filt[,x]), mc.cores=num_cores)
betas.filt.norm <- do.call(cbind, betas.filt.norm_list)  # A matrix (not RatioSet)
rm(mSet.filt, betas.filt.norm_list)
print("BMIQ normalization done. Memory used:")
print(pryr::mem_used())
p6 <- proc.time()
print(p6-p5)
save("betas.filt.norm", file=paste0("../int/betas.filt.norm_", job_idx, ".RData"))  # Save as .RData with associated index no.

print("BMIQ-normalized matrix saved.")
print(pryr::mem_used())
p7 <- proc.time()
print(p7-p6)

# cellCounts <- estimateCellCounts(rgSet)
WBC_relevant_sites <- betas.filt.norm[rownames(betas.filt.norm) %in% 
                                              rownames(FlowSorted.Blood.450k.JaffeModelPars),]
cellCounts <- projectMix(WBC_relevant_sites, FlowSorted.Blood.450k.JaffeModelPars)
save("cellCounts", file=paste0("../int/cellCounts_", job_idx, ".RData"))  # Save as .RData with associated index no.
print("Cell counts matrix saved.")
print(pryr::mem_used())
p8 <- proc.time()
print(p8-p7)
