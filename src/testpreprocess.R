# # system.time({
# #   idx_split <- splitIndices(50, 16)
# #   rgList <- parLapply(cl, idx_split, function(idxs) read.metharray.exp(targets=targets[idxs,]))
# #   greens <- do.call(cbind, lapply(rgList, getGreen))
# #   reds <- do.call(cbind, lapply(rgList, getRed))
# #   fullRGSet <- RGChannelSet(Green=greens, Red=reds, annotation=annotation(rgList[[1]]))
# # })
# # pryr::mem_used()
# # mSet <- preprocessRaw(fullRGset)
# # pryr::mem_used()
# # bmiq_list <- 
# # ---------------
# 
suppressMessages(silent <- lapply(c("tidyverse","minfi","wateRmelon","RefFreeEWAS","FlowSorted.Blood.450k"),
                                    library, character.only=T))

# ## Sample sheet
# targets <- read.metharray.sheet(base="../data/", pattern="^sample_sheet4.csv$") %>%
#   dplyr::filter(Basename!="character(0)",  # Seems to appear in read.metharray.sheet output when a channel is missing
#                 Sample_Name!=21084)  # This sample is missing red channel .idat
# targets <- targets[1:700,]
# 
args <- commandArgs(trailingOnly=T)
num_cores <- as.numeric(args[1])
# 
# ## Set up group of workers for parallelization
# p1 <- proc.time()
# cl <- makeCluster(num_cores)
# clusterExport(cl, "targets")
# silent <- clusterEvalQ(cl, library(minfi))
# rgSet_list <- parLapply(cl, splitIndices(nrow(targets), num_cores),  # Splits samples into equally-sized chunks
#                         function(idx_set) read.metharray.exp(targets=targets[idx_set,]))
# # stopCluster(cl)
# greens <- do.call(cbind, lapply(rgSet_list, getGreen))
# reds <- do.call(cbind, lapply(rgSet_list, getRed))
# rgSet <- RGChannelSet(Green=greens, Red=reds, annotation=annotation(rgSet_list[[1]]))
# rm(rgSet_list, greens, reds)
# print("RGChannelSet constructed.")
# p2 <- proc.time()
# print(p2-p1)
# pryr::mem_used()
# save(rgSet, file="../int/toyrgSet.RData")
# # ## Detection-based QC
# # detP <- detectionP(rgSet)
# # print("Detection-based QC run. Mem used so far:")
# # pryr::mem_used()
# # p4 <- proc.time()
# # print(p4-p3)
# # keep_samples <- colSums(detP>0.01)/nrow(detP) < 0.1  # Don't want >10% of probes whose detection failed at p=0.01
# # detP <- detP[,keep_samples]  # Don't look at probes until bad samples have been removed
# # bad_probes <- rownames(detP)[rowSums(detP>0.01)/ncol(detP) > 0.01]  # Save these and filter out later
# # save("bad_probes", file=paste0("../int/bad_probes_", job_idx, ".RData"))
# # mSet.filt <- mSet[,keep_samples]  # Remove low-quality samples
# # print(paste("Removed", sum(!keep_samples), "samples"))
# # mSet.filt <- dropMethylationLoci(mSet.filt)  # Drop SNP probes and CpH probes
# # rm(detP, mSet)
# # print("QC complete. Mem used so far:")
# # pryr::mem_used()
# # p5 <- proc.time()
# # print(p5-p4)

load("../int/toyrgSet.RData")
# ## foreach method
# cl <- makeCluster(num_cores)
# library(doParallel)
# registerDoParallel(cl)
# mSet_list <- foreach(i=splitIndices(ncol(rgSet), num_cores)) %dopar%
#   function(idx_set) preprocessNoob(rgSet[,idx_set])
# print("mSet constructed using preprocessNoob -- didn't die from memory!")

## parLapply method
cl <- makeCluster(16, useXDR=F)
clusterExport(cl, "rgSet")
silent <- clusterEvalQ(cl, library(minfi))
# silent <- clusterEvalQ(cl, gc())
# splitrgSet <- lapply(clusterSplit(cl, 1:ncol(rgSet)), function(idx) rgSet[,idx])
# rm(rgSet)
mSet_list <- parLapply(cl, splitIndices(ncol(rgSet), num_cores),
                       function(idx) preprocessNoob(rgSet[,idx], dyeMethod="single"))
stopCluster(cl)
print("mSet constructed using preprocessNoob -- didn't die from memory!")


# cellCounts <- estimateCellCounts(rgSet)
WBC_relevant_sites <- betas.filt.norm[rownames(betas.filt.norm) %in%
                                        rownames(FlowSorted.Blood.450k.JaffeModelPars),]
cellCounts <- projectMix(WBC_relevant_sites, FlowSorted.Blood.450k.JaffeModelPars)
# save("cellCounts", file=paste0("../int/cellCounts_", job_idx, ".RData"))  # Save as .RData with associated index no.
# print("Cell counts matrix saved.")
# print(pryr::mem_used())
# p8 <- proc.time()
# print(p8-p7)




