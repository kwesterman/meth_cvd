suppressMessages(lapply(c("tidyverse","minfi"), library, character.only=T))

args <- commandArgs(trailingOnly=T)
num_cores <- as.numeric(args[1])-1  # Get number of available cores as script argument

## Read in preprocessed beta values
# # For toy:
# betas.filt.norms <- lapply(1:7, function(x) betas.filt.norm[,(2*x-1):(2*x)])
# betas.filt.norm <- do.call(cbind, lapply(betas.filt.norms, getBeta))

# For cluster
betas.filt.norm_files <- grep("betas.filt.norm_[0-9]+.RData", list.files("../int/"), value=T)
betas.filt.norms <- mclapply(paste0("../int/", betas.filt.norm_files), 
                             function(x) {load(x); betas.filt.norm}, 
                             mc.cores=num_cores, mc.preschedule=F)
betas.filt.norm <- do.call(cbind, betas.filt.norms)

## Retrieve technical sample covariates
sampleData <- read.metharray.sheet(base="../data/", pattern="^sample_sheet4.csv$") %>%
  dplyr::mutate(key=paste(Slide, Array, sep="_")) %>%
  dplyr::slice(match(colnames(betas.filt.norm), key)) %>% 
  dplyr::mutate(row=factor(substring(Array, 3, 3)), 
                col=factor(substring(Array, 6, 6)),
                Slide=factor(Slide)) %>%
  dplyr::select(key, Slide, row, col)
stopifnot(nrow(sampleData) == ncol(betas.filt.norm))

## QC Plots (data too big to be efficient/useful? useless after normalization?)
# jpeg("../output/densityPlot.jpg")
# densityPlot(betas.filt.norm, sampGroups=NULL)  # Subsitute label vector for any sample group of interest here
# dev.off()
# jpeg("../output/densityBeanPlot.jpg")
# densityBeanPlot(betas.filt.norm, sampGroups=NULL)  # Subsitute label vector for any sample group of interest here
# dev.off()

## Filter out bad probes as identified earlier
badProbe_files <- grep("bad_probes_[0-9]+.RData", list.files("../int/"), value=T)
print(badProbe_files)
badProbes_list <- mclapply(paste0("../int/", badProbe_files), 
                            function(x) {load(x); bad_probes}, mc.cores=num_cores)
badProbes <- unique(unlist(badProbes_list))
betas.filt.norm <- betas.filt.norm[!(rownames(betas.filt.norm) %in% badProbes),]
print("Probes filtered.")

## Generate M-values (= logit2(betas)) for subsequent analysis
Mvals <- logit2(betas.filt.norm)
print("Mvals generated")
rm(betas.filt.norm_files, betas.filt.norms, betas.filt.norm)
save("Mvals", file="../int/Mvals.RData")
print("M-values generated and saved.")

## Run CPA (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
rgSet_files <- grep("rgSet_[0-9]+.RData", list.files("../int/"), value=T)
source("Lehne2015.R")
CP_mat_list <- lapply(paste0("../int/", rgSet_files), function(x) {load(x); get_CPs(rgSet)})
# cl <- makeCluster(num_cores)
# clusterExport(cl, "CPA")
# clusterEvalQ(cl, library(minfi))
# CP_mat_list <- parLapply(cl, paste0("../int/", rgSet_files), function(x) {load(x); get_CPs(rgSet)})
# stopCluster(cl)
# CP_mat_list <- mclapply(,  # Use adapted code from Lehne 2015 to produce 235 x n_samples matrix of control probe intensities
#                    , mc.cores=num_cores)
CP_mat <- do.call(cbind, CP_mat_list)
print("Control probe matrix constructed. Performing CPA...")
print(dim(t(CP_mat)))
pca <- prcomp(na.omit(t(CP_mat)))
CP_PCs <- pca$x
colnames(CP_PCs) = paste(colnames(CP_PCs), '_cp', sep='')
print("CPA done.")

save("CP_PCs", file="../int/CP_PCs.RData")
print("CPA performed and results saved.")

# ## ComBat batch correction on M-values: avoid for the moment
# # Estimated ~15 secs per sample
# betas.filt.norm_Mvals <- logit2(betas.filt.norm)
# system.time(a <- ComBat(dat=betas.filt.norm_Mvals, batch=targets$Slide))


# technical_regression <- function(i, sampleData, betas.filt.norm_Mvals) {
#   sampleData$Mvals <- betas.filt.norm_Mvals[i,]
#   tryCatch(lm(Mvals~Slide+row+col, data=sampleData)$residuals,
#            error = function(e) rep(NA, nrow(sampleData)))
# }
# print(pryr::mem_used())
# print("starting regressions")
# num_regs <- nrow(betas.filt.norm_Mvals)
# print(num_regs)
# cl <- makeCluster(num_cores)
# resids <- parLapply(cl, setNames(1:num_regs, rownames(betas.filt.norm_Mvals)[1:num_regs]), 
#                     technical_regression, sampleData, betas.filt.norm_Mvals)
# stopCluster(cl)
# # resids <- mclapply(setNames(1:num_regs, rownames(betas.filt.norm)[1:num_regs]), technical_regression, mc.cores=8)
# # resids <- lapply(setNames(1:num_regs, rownames(betas.filt.norm)[1:num_regs]), technical_regression)
# print("regressions complete")
# pryr::mem_used()
# meth_adj_Mvals <- do.call(rbind, resids)
# pryr::mem_used()
# print(sort(sapply(ls(), function(x) object.size(get(x)))))
# rm(resids)
# colnames(meth_adj_Mvals) <- colnames(betas.filt.norm_Mvals)
# print("saving")
# save("meth_adj_Mvals", file="../int/meth_adj.RData")
