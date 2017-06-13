
lapply(c("knitr","tidyverse","minfi","wateRmelon","RefFreeEWAS","FlowSorted.Blood.450k"),
       library, character.only=T)

dataDir <- "../data/"  # At the moment, this only contains a few samples for development
targets <- read.metharray.sheet(base=dataDir, pattern="^sample_sheet4.csv$") %>%
  dplyr::filter(Basename!="character(0)",
                Sample_Name!=21084)  # This sample is missing .idat for red channel
rgSet <- read.metharray.exp(targets=targets[1:5,], extended=T, verbose=T)

mSet.filt <- pfilter(mn=rgSet, pnthresh=0.05, perc=5, pthresh=5, useNoob=F)  # Samples and probes must have <5% of representatives with detection P>0.05

my_filter <- function(rgSet, useNoob=F, pval_thresh=0.05, probe_thresh=0.01, sample_thresh=0.01) {
  detP <- detectionP(rgSet)  # results in n_features x n_samples matrix of detection p-values
  keep_probes <- rowSums(detP>pval_thresh)/ncol(detP) < probe_thresh
  keep_samples <- colSums(detP>pval_thresh)/nrow(detP) < sample_thresh
  preprocessed <- {if (useNoob) preprocessNoob(rgSet) else preprocessRaw(rgSet)}
  print(paste("Removed", sum(!keep_samples), "samples"))
  print(paste("Removed", sum(!keep_probes), "probes"))
  stopifnot(dim(detP)==dim(preprocessed))
  preprocessed[keep_probes, keep_samples]
}
my_filter(nonext)

mSet.filt.bmiq_matrix <- BMIQ(mSet.filt)
mSet.filt.bmiq <- RatioSet(Beta=mSet.filt.bmiq_matrix)

data(FlowSorted.Blood.450k.compTable)
cellRef <- as.matrix(FlowSorted.Blood.450k.compTable)[,c('CD8T','CD4T','NK','Bcell','Mono','Gran')]  # Only needs cols representing actual cell types
commonCpGs <- intersect(rownames(mSet.filt.bmiq), rownames(cellRef))  # Which CpG sites are common to my set and the cell type reference?
estCellCounts <- projectMix(Y=mSet.filt.bmiq_matrix[commonCpGs,], Xmat=cellRef[commonCpGs,])

targets_withKey <- mutate(targets, Sample_pKey=paste(Slide, Sample_Well, sep="_"))
estCellCounts_df <- rownames_to_column(as.data.frame(estCellCounts), var="Sample_pKey")
meth_covData <- inner_join(targets_withKey, estCellCounts_df, by="Sample_pKey")

save("mSet.filt.bmiq", file="../data/preprocessed_methylation_data.Rdata")
save("meth_covData", file="../data/methylation_covariates.Rdata")


