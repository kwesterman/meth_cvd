suppressMessages(lapply(c("tidyverse","minfi","wateRmelon","RefFreeEWAS","FlowSorted.Blood.450k"),
       library, character.only=T))

print("Before loading anything:")
load("../data/rgSet_full.RData")  # Loads object "rgSet"
print(paste("rgSet loaded. Memory used: ", pryr::mem_used()))

my_filter <- function(rgSet, useNoob=F, pval_thresh=0.05, probe_thresh=0.01, sample_thresh=0.01) {
  ## Perform QC (sample and probe filtering) and optional background correction
  ## Input: RGChannelSet
  ## Output: preprocessed MethylSet
  ## Noob option performs background correction and dye bias adjustment
  detP <- detectionP(rgSet)  # results in n_features x n_samples matrix of detection p-values
  keep_probes <- rowSums(detP>pval_thresh)/ncol(detP) < probe_thresh
  keep_samples <- colSums(detP>pval_thresh)/nrow(detP) < sample_thresh
  preprocessed <- {if (useNoob) preprocessNoob(rgSet) else preprocessRaw(rgSet)}
  print(paste("Removed", sum(!keep_samples), "samples"))
  print(paste("Removed", sum(!keep_probes), "probes"))
  stopifnot(dim(detP)==dim(preprocessed))
  preprocessed[keep_probes, keep_samples]
}

p1 <- proc.time()

mSet.filt <- my_filter(rgSet, useNoob=F)
print("Filtering done.")
print(pryr::mem_used())
p2 <- proc.time()
print(p2-p1)


mSet.filt.bmiq_matrix <- BMIQ(mSet.filt)
mSet.filt.bmiq <- RatioSet(Beta=mSet.filt.bmiq_matrix)
print("BMIQ done.")
print(pryr::mem_used())
p3 <- proc.time()
print(p3-p2)

data(FlowSorted.Blood.450k.compTable)
cellRef <- as.matrix(FlowSorted.Blood.450k.compTable)[,c('CD8T','CD4T','NK','Bcell','Mono','Gran')]  # Only needs cols representing actual cell types
commonCpGs <- intersect(rownames(mSet.filt.bmiq), rownames(cellRef))  # Which CpG sites are common to my set and the cell type reference?
estCellCounts <- projectMix(Y=mSet.filt.bmiq_matrix[commonCpGs,], Xmat=cellRef[commonCpGs,])
load("../data/targets.RData")
targets_withKey <- mutate(targets, Sample_pKey=paste(Slide, Sample_Well, sep="_"))
estCellCounts_df <- rownames_to_column(as.data.frame(estCellCounts), var="Sample_pKey")
meth_covData <- inner_join(targets_withKey, estCellCounts_df, by="Sample_pKey")
print("covData matrix creation done.")
print(pryr::mem_used())
p4 <- proc.time()
print(p4-p3)

save("mSet.filt.bmiq", file="../data/preprocessed_methylation_data.Rdata")
print("BMIQ-normalized matrix saved.")
print(pryr::mem_used())
p5 <- proc.time()
print(p5-p4)
save("meth_covData", file="../data/methylation_covariates.Rdata")