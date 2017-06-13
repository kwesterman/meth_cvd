suppressMessages(silent <- lapply(c("tidyverse","minfi","wateRmelon","doParallel"), library, character.only=T))

num_cores <- min(detectCores(), 8)
print(paste("Number of cores:", num_cores))

load("../int/rgSet.RData")  # Load RGChannelSet object from previous step in pipeline

## Noob procedure for background correction and dye bias adjustment
print("Preprocessing using Noob method...")

rgSet.split <- lapply(splitIndices(ncol(rgSet),32), function(idx_set) rgSet[,idx_set]) 
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
mSet_list <- foreach(rgChunk=rgSet.split, .packages="minfi") %dopar% preprocessNoob(rgChunk)
stopCluster(cl)
mSet <- MethylSet(Meth=do.call(cbind, lapply(mSet_list, getMeth)),  # Construct MethylSet from mSet_list chunks
                  Unmeth=do.call(cbind, lapply(mSet_list, getUnmeth)),
                  phenoData=do.call(Biobase::combine, lapply(mSet_list, phenoData)))
rm(rgSet, rgSet.split, mSet_list); invisible(gc())


## QC/Filtering
print("QC/Filtering steps...")
load("../int/detP.RData")  # Pre-calculated detection p-values
keepSamples <- colSums(detP>0.01)/nrow(detP) < 0.1  # Don't want >10% of probes whose detection failed at p=0.01
undetectedProbes <- rowSums(detP>0.01)/ncol(detP) < 0.1  # Don't want >10% of samples whose detection failed at p=0.01
ann450k <- getAnnotation(mSet)
sexChromProbes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")]  # Probes in sex chromosomes
crossReactiveProbesDF <- read.csv("../data/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=F)  # From Chen 2013
snpProbes <- grep("^rs", featureNames(mSet), value=T)  # Probes measuring SNPs (may be none)
chProbes <- grep("^ch\\.", featureNames(mSet), value=T)  # Probes measuring CpH (non-CpG) methylation
snpInfo <- getSnpInfo(mSet)
probesWithSNPs <- rownames(snpInfo)[!is.na(snpInfo$CpG_maf) | !is.na(snpInfo$SBE_maf)]  # Probes with SNPs at CpG site or single-base extension site
keepProbes <- !(rownames(mSet) %in% c(undetectedProbes, sexChromProbes, crossReactiveProbesDF$TargetID,
                                      snpProbes, chProbes, probesWithSNPs))
mSet.filt <- mSet[keepProbes,keepSamples]
print(paste("QC: Removed", sum(!keepSamples), "samples."))
print(paste("Dimensions of filtered MethylSet:", dim(mSet.filt)[1], "x", dim(mSet.filt)[2]))
jpeg("../output/qcPlot.jpg"); plotQC(getQC(mSet.filt)); dev.off()  # Plot of Unmeth vs. Meth median intensities
rm(detP, mSet); invisible(gc())


## Within-sample normalization using BMIQ
print("Normalization using BMIQ...")

mSet.filt.split <- lapply(splitIndices(ncol(mSet.filt),32), function(idx_set) mSet.filt[,idx_set]) 
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
betas.filt.norm <- foreach(mChunk=mSet.filt.split, .combine=cbind,.packages=c("minfi","wateRmelon")) %dopar% 
  BMIQ(mChunk)
stopCluster(cl)
print("BMIQ normalization done.")
rm(mSet.filt, mSet.filt.split); invisible(gc())

print("Saving M-value matrix...")
Mvals <- logit2(betas.filt.norm)
save("Mvals", file="../int/Mvals.RData")
