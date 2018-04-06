suppressMessages(silent <- lapply(c("tidyverse","minfi","wateRmelon","doParallel"), library, character.only=T))

datasets <- commandArgs(TRUE)  # Create list of datasets to use from command-line arguments

numCores <- min(detectCores(), 16)
qc_thresh <- list(whi=c(mMed=11, uMed=10.5),
                  fhs=c(mMed=10, uMed=10))

run_noob <- function(rgSet, numCores) {
  # Noob procedure for background correction and dye bias adjustment
  print("...preprocessing using Noob method...")
  
  rgSet.split <- lapply(splitIndices(ncol(rgSet),2*numCores), function(idx_set) rgSet[,idx_set])
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  mSet.list <- foreach(rgChunk=rgSet.split, .packages="minfi") %dopar% 
    preprocessNoob(rgChunk)
  mSet <- MethylSet(Meth=do.call(cbind, lapply(mSet.list, getMeth)),
                    Unmeth=do.call(cbind, lapply(mSet.list, getUnmeth)),
                    phenoData=do.call(combine, lapply(mSet.list, phenoData)))
  stopCluster(cl)
  
  mSet
}

run_sample_qc <- function(dataset, mSet) {
  # Sample QC
  print("...sample QC...")
  detP <- readRDS(paste0("../int/detP_", dataset, ".rds"))
  lowDetection <- colSums(detP>1e-16)/nrow(detP) > 0.1  # Don't want >10% of probes whose detection failed at p=0.01
  medIntensities <- getQC(mSet)
  lowIntensity <- (medIntensities$mMed < qc_thresh[[dataset]]["mMed"] | 
                     medIntensities$uMed < qc_thresh[[dataset]]["uMed"])  # Informed by visual inspection of U vs. M intensity plot
  sexMismatch <- pData(mSet)$sex!=getSex(mapToGenome(mSet))$predictedSex
  keepSamples <- !(lowDetection | lowIntensity | sexMismatch)
  mSet.qc <- mSet[,keepSamples]
  print(paste0("...QC: Removed ", sum(!keepSamples), " samples (", sum(lowDetection), " detection, ", 
              sum(lowIntensity), " intensity, ", sum(sexMismatch), " sex mismatch)."))
  
  mSet.qc
}

make_dplots <- function(betas, types, ttl) {
  # Type I/II comparison plot
  plot(0, xlim=c(-0.2,1.2), ylim=c(0,10), xlab="Beta values", ylab="Density", main=ttl)
  typeI_densities <- apply(betas[types=="I",], 2, density)
  for (samp in typeI_densities) lines(samp, col="red")
  typeII_densities <- apply(betas[types=="II",], 2, density)
  for (samp in typeII_densities) lines(samp, col="green")
  legend("topright", c("Type I", "Type II"), col=c("red","green"), lty=1)
}

run_bmiq <- function(dataset, mSet.qc) {
  # Within-sample normalization using BMIQ
  print("...normalization using BMIQ...")
  
  random1k <- sample.int(ncol(mSet.qc), size=1000)
  jpeg(paste0("../output/probeTypePlot_prenormalization_", dataset, ".jpg"))
  make_dplots(getBeta(mSet.qc[,random1k]), getAnnotation(mSet.qc)$Type, ttl="Pre-normalization")
  dev.off()
  
  mSet.qc.split <- lapply(splitIndices(ncol(mSet.qc),2*numCores), function(idx_set) mSet.qc[,idx_set])
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  betas.qc.norm_list <- foreach(mChunk=mSet.qc.split, .packages=c("minfi","wateRmelon")) %dopar% 
    BMIQ(mChunk)
  betas.qc.norm <- do.call(cbind, betas.qc.norm_list)
  stopCluster(cl)
  
  jpeg(paste0("../output/probeTypePlot_postnormalization_", dataset, ".jpg"))
  make_dplots(betas.qc.norm[,random1k], getAnnotation(mSet.qc)$Type, ttl="Post-normalization")
  dev.off()
  
  print("BMIQ normalization done.")
  betas.qc.norm
}

filter_probes <- function(mSet.qc, betas.qc.norm) {
  # Probe filtering
  print("...probe filtering...")
  detP <- readRDS(paste0("../int/detP_", dataset, ".rds"))
  undetectedProbes <- rownames(detP)[rowSums(detP>1e-16)/ncol(detP) > 0.1]  # Don't want >10% of samples whose detection failed at p=0.01
  ann450k <- getAnnotation(mSet.qc)
  sexChromProbes <- ann450k$Name[ann450k$chr %in% c("chrX","chrY")]  # Probes in sex chromosomes
  crossReactiveProbesDF <- read.csv("../data/literature/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=F)  # From Chen 2013
  snpProbes <- grep("^rs", featureNames(mSet.qc), value=T)  # Probes measuring SNPs (may be none)
  chProbes <- grep("^ch\\.", featureNames(mSet.qc), value=T)  # Probes measuring CpH (non-CpG) methylation
  snpInfo <- getSnpInfo(mSet.qc)
  probesWithSNPs <- rownames(snpInfo)[!is.na(snpInfo$CpG_maf) | !is.na(snpInfo$SBE_maf)]  # Probes with SNPs at CpG site or single-base extension site
  keepProbes <- !(rownames(mSet.qc) %in% c(undetectedProbes, sexChromProbes, crossReactiveProbesDF$TargetID,
                                           snpProbes, chProbes, probesWithSNPs))
  betas.qc.norm.filt <- betas.qc.norm[keepProbes,]
  print(paste("...filtering: removed", sum(!keepProbes), "probes."))
  
  betas.qc.norm.filt
}

for (dataset in datasets) {
  print(paste0("Dataset: ", dataset))
  rgSet <- readRDS(paste0("../int/rgSet_", dataset, ".rds"))  # Load RGChannelSet object from previous step in pipeline
  mSet <- run_noob(rgSet, numCores)
  mSet.qc <- run_sample_qc(dataset, mSet)
  betas.qc.norm <- run_bmiq(dataset, mSet.qc)
  betas.qc.norm.filt <- filter_probes(mSet.qc, betas.qc.norm)
  print(paste("...dimensions of final methylation set:", paste(dim(betas.qc.norm.filt), collapse=" x ")))
  print("...saving methylation matrix...")
  saveRDS(betas.qc.norm.filt, file=paste0("../int/betas.qc.norm.filt_", dataset, ".rds"), compress=F)
  rm(rgSet, mSet, mSet.qc, betas.qc.norm, betas.qc.norm.filt)
}

