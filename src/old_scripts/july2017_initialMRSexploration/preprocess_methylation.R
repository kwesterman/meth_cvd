suppressMessages(silent <- lapply(c("tidyverse","minfi","wateRmelon","doParallel"), library, character.only=T))

num_cores <- min(detectCores(), 8)
print(paste("Number of cores:", num_cores))

load("../int/rgSet.RData")  # Load RGChannelSet object from previous step in pipeline

## Noob procedure for background correction and dye bias adjustment
print("Preprocessing using Noob method...")

rgSet.split <- lapply(splitIndices(ncol(rgSet),32), function(idx_set) rgSet[,idx_set]) 
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
mSet <- foreach(rgChunk=rgSet.split, .packages="minfi", .combine="combine") %dopar% 
  preprocessNoob(rgChunk)
stopCluster(cl)
# mSet <- MethylSet(Meth=do.call(cbind, lapply(mSet_list, getMeth)),  # Construct MethylSet from mSet_list chunks
#                   Unmeth=do.call(cbind, lapply(mSet_list, getUnmeth)),
#                   phenoData=do.call(Biobase::combine, lapply(mSet_list, phenoData)))
rm(rgSet, rgSet.split); invisible(gc())

## Sample QC
print("Sample QC...")
load("../int/detP.RData")
lowDetection <- colSums(detP>0.01)/nrow(detP) > 0.1  # Don't want >10% of probes whose detection failed at p=0.01
medIntensities <- getQC(mSet)
lowIntensity <- medIntensities$mMed < 10 | medIntensities$uMed < 10  # Informed by visual inspection of U vs. M intensity plot
sexMismatch <- pData(mSet)$sex!=getSex(mapToGenome(mSet))$predictedSex
keepSamples <- !(lowDetection | lowIntensity | sexMismatch)
mSet.qc <- mSet[,keepSamples]
print(paste("QC: Removed", sum(!keepSamples), "samples."))
print(paste(sum(lowDetection), "detection,", sum(lowIntensity), "intensity,", sum(sexMismatch), "sex mismatch"))
rm(mSet); invisible(gc())

## Type I/II comparison plot (pre-normalization)
make_dplots <- function(betas, types, ttl) {
  plot(0, xlim=c(-0.2,1.2), ylim=c(0,10), xlab="Beta values", ylab="Density", main=ttl)
  typeI_densities <- apply(betas[types=="I",], 2, density)
  for (samp in typeI_densities) lines(samp, col="red")
  typeII_densities <- apply(betas[types=="II",], 2, density)
  for (samp in typeII_densities) lines(samp, col="green")
  legend("topright", c("Type I", "Type II"), col=c("red","green"), lty=1)
}
jpeg("../output/probeTypePlot_prenormalization.jpg")
make_dplots(getBeta(mSet.qc), getAnnotation(mSet.qc)$Type, ttl="Pre-normalization")
dev.off()

## Within-sample normalization using BMIQ
print("Normalization using BMIQ...")

mSet.qc.split <- lapply(splitIndices(ncol(mSet.qc),32), function(idx_set) mSet.qc[,idx_set]) 
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
betas.qc.norm <- foreach(mChunk=mSet.qc.split, .combine=cbind,.packages=c("minfi","wateRmelon")) %dopar% 
  BMIQ(mChunk)
stopCluster(cl)
print("BMIQ normalization done.")
rm(mSet.qc.split); invisible(gc())

## Type I/II comparison plot (post-normalization)
jpeg("../output/probeTypePlot_postnormalization.jpg")
make_dplots(betas.qc.norm, getAnnotation(mSet.qc)$Type, ttl="Post-normalization")
dev.off()

## Probe filtering
print("Probe filtering...")
# undetectedProbes <- rowSums(detP>0.01)/ncol(detP) < 0.1  # Don't want >10% of samples whose detection failed at p=0.01
ann450k <- getAnnotation(mSet.qc)
sexChromProbes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")]  # Probes in sex chromosomes
crossReactiveProbesDF <- read.csv("../data/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=F)  # From Chen 2013
snpProbes <- grep("^rs", featureNames(mSet.qc), value=T)  # Probes measuring SNPs (may be none)
chProbes <- grep("^ch\\.", featureNames(mSet.qc), value=T)  # Probes measuring CpH (non-CpG) methylation
snpInfo <- getSnpInfo(mSet.qc)
probesWithSNPs <- rownames(snpInfo)[!is.na(snpInfo$CpG_maf) | !is.na(snpInfo$SBE_maf)]  # Probes with SNPs at CpG site or single-base extension site
keepProbes <- !(rownames(mSet.qc) %in% c(sexChromProbes, crossReactiveProbesDF$TargetID,
                                      snpProbes, chProbes, probesWithSNPs))
betas.qc.norm.filt <- betas.qc.norm[keepProbes,]
print(paste("Filtering: Removed", sum(!keepProbes), "probes."))
print(paste("Dimensions of final methylation set:", dim(betas.qc.norm.filt)[1], "x", dim(betas.qc.norm.filt)[2]))
rm(detP); invisible(gc())


print("Saving M-value matrix...")
Mvals <- logit2(betas.qc.norm.filt)
save("Mvals", file="../int/Mvals.RData")
