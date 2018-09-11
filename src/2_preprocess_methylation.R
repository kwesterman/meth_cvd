suppressMessages(silent <- lapply(
  c("tidyverse", "minfi", "wateRmelon", "doParallel"), 
  library, character.only=T))

datasets <- commandArgs(TRUE)  # List of datasets to use from arguments

num_cores <- min(detectCores(), 16)
qc_thresh <- list(whi=c(m_med=11, u_med=10.5),
                  fhs=c(m_med=10, u_med=10))

run_noob <- function(rg_set, num_cores) {
  # Noob procedure for background correction and dye bias adjustment
  print("...preprocessing using Noob method...")
  
  rg_set_split <- lapply(splitIndices(ncol(rg_set), 2*num_cores), 
                         function(idx_set) rg_set[, idx_set])
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  m_set.list <- foreach(rgChunk=rg_set_split, .packages="minfi") %dopar% 
    preprocessNoob(rgChunk)
  m_set <- MethylSet(Meth=do.call(cbind, lapply(m_set.list, getMeth)),
                    Unmeth=do.call(cbind, lapply(m_set.list, getUnmeth)),
                    phenoData=do.call(combine, lapply(m_set.list, phenoData)))
  stopCluster(cl)
  
  m_set
}

run_sample_qc <- function(dataset, m_set) {
  # Sample QC
  print("...sample QC...")
  detP <- readRDS(paste0("../int/detP_", dataset, ".rds"))
  low_detection <- colSums(detP>1e-16) / nrow(detP) > 0.1  # Don't want >10% of probes whose detection failed at p=1e-16
  med_intensities <- getQC(m_set)
  low_intensity <- (med_intensities$mMed < qc_thresh[[dataset]]["mMed"] | 
                      med_intensities$uMed < qc_thresh[[dataset]]["uMed"])  # Informed by visual inspection of U vs. M intensity plot
  sex_mismatch <- pData(m_set)$sex != getSex(mapToGenome(m_set))$predictedSex
  keep_samples <- !(low_detection | low_intensity | sex_mismatch)
  m_set_qc <- m_set[, keep_samples]
  print(paste0("...QC: Removed ", sum(!keep_samples), " samples (",
               sum(low_detection), " detection, ", 
               sum(low_intensity), " intensity, ", 
               sum(sex_mismatch), " sex mismatch)."))
  m_set_qc
}

make_dplots <- function(betas, types, ttl) {
  # Type I/II comparison plot
  plot(0, xlim=c(-0.2, 1.2), ylim=c(0, 10), 
       xlab="Beta values", ylab="Density", main=ttl)
  typeI_densities <- apply(betas[types == "I", ], 2, density)
  for (samp in typeI_densities) lines(samp, col="red")
  typeII_densities <- apply(betas[types == "II", ], 2, density)
  for (samp in typeII_densities) lines(samp, col="green")
  legend("topright", c("Type I", "Type II"), col=c("red", "green"), lty=1)
}

run_bmiq <- function(dataset, m_set_qc) {
  # Within-sample normalization using BMIQ
  print("...normalization using BMIQ...")
  
  random_1k <- sample.int(ncol(m_set_qc), size=1000)
  jpeg(paste0("../output/probeTypePlot_prenormalization_", dataset, ".jpg"))
  make_dplots(getBeta(m_set_qc[, random_1k]), getAnnotation(m_set_qc)$Type, 
              ttl="Pre-normalization")
  dev.off()
  
  m_set_qc_split <- lapply(splitIndices(ncol(m_set_qc), 2 * num_cores), 
                           function(idx_set) m_set_qc[, idx_set])
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  betas_qc_norm_list <- foreach(m_chunk=m_set_qc_split, 
                                .packages=c("minfi","wateRmelon")) %dopar% 
    BMIQ(m_chunk)
  betas_qc_norm <- do.call(cbind, betas_qc_norm_list)
  stopCluster(cl)
  
  jpeg(paste0("../output/probeTypePlot_postnormalization_", dataset, ".jpg"))
  make_dplots(betas_qc_norm[, random_1k], getAnnotation(m_set_qc)$Type, 
              ttl="Post-normalization")
  dev.off()
  
  print("BMIQ normalization done.")
  betas_qc_norm
}

filter_probes <- function(m_set_qc, betas_qc_norm) {
  # Probe filtering
  print("...probe filtering...")
  detP <- readRDS(paste0("../int/detP_", dataset, ".rds"))
  undetected_probes <- rownames(detP)[rowSums(detP>1e-16) / ncol(detP) > 0.1]  # Don't want >10% of samples whose detection failed at p=1e-16
  anno_450k <- getAnnotation(m_set_qc)
  sex_chrom_probes <- anno_450k$Name[anno_450k0k$chr %in% c("chrX", "chrY")]  # Probes in sex chromosomes
  cross_reactive_probes_df <- read.csv(
    "../data/literature/48639-non-specific-probes-Illumina450k.csv", 
    stringsAsFactors=F)  # From Chen 2013
  snp_probes <- grep("^rs", featureNames(m_set_qc), value=T)  # Probes measuring SNPs (may be none)
  ch_probes <- grep("^ch\\.", featureNames(m_set_qc), value=T)  # Probes measuring CpH (non-CpG) methylation
  snp_info <- getSnpInfo(m_set_qc)
  probes_with_SNPs <- rownames(snp_info)[(
    !is.na(snp_info$CpG_maf) | !is.na(snp_info$SBE_maf))]  # Probes with SNP at CpG or single-base extension site
  keep_probes <- !(rownames(m_set_qc) %in% c(
    undetected_probes, sex_chrom_probes, cross_reactive_probes_df$TargetID,
    snp_probes, ch_probes, probes_with_SNPs))
  betas_qc_norm_filt <- betas_qc_norm[keep_probes,]
  print(paste("...filtering: removed", sum(!keep_probes), "probes."))
  
  betas_qc_norm_filt
}

for (dataset in datasets) {
  print(paste0("Dataset: ", dataset))
  rg_set <- readRDS(paste0("../int/rg_set_", dataset, ".rds"))  # Load RGChannelSet object from previous step in pipeline
  m_set <- run_noob(rg_set, num_cores)
  m_set_qc <- run_sample_qc(dataset, m_set)
  betas_qc_norm <- run_bmiq(dataset, m_set_qc)
  betas_qc_norm_filt <- filter_probes(m_set_qc, betas_qc_norm)
  print(paste("...dimensions of final methylation set:", 
              paste(dim(betas_qc_norm_filt), collapse=" x ")))
  print("...saving methylation matrix...")
  saveRDS(betas_qc_norm_filt, 
          file=paste0("../int/betas_qc_norm_filt_", dataset, ".rds"), 
          compress=F)
  rm(rg_set, m_set, m_set_qc, betas_qc_norm, betas_qc_norm_filt)
}

