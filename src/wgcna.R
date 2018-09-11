library(tidyverse)
library(WGCNA)
library(survival)

betas_whi <- readRDS("../int/betas.qc.norm.filt_whi.rds")
beta_variances <- apply(betas_whi, 1, var)

enableWGCNAThreads(8)  # use enableWGCNAThreads() instead of allow...()


# Subset of individuals for threshold choice and pre-clustering
set.seed(123)
sample_subset <- sample(1:ncol(betas_whi), 100)  
meth_dat_subset <- t(betas_whi[beta_variances > quantile(beta_variances, 0.5),
                               sample_subset])

# Analysis for power choice
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
system.time(sft <- pickSoftThreshold(meth_dat_subset, 
                                     powerVector=powers, verbose=5))
## soft threshold choice on 100 samples for top 50% variable genes took 3.5 hrs
saveRDS(sft, file="../int/wgcna/whi_soft_threshold.rds")

# Code from tutorial to assess proper soft threshold choice
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main=paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels=powers, cex=cex1, col="red")
## --> based on visual analysis, preliminary choice of beta=6

# Pre-clustering using projective K-means
system.time(pkm <- projectiveKMeans(t(betas_whi[, sample_subset])))
## pkm on 100 samples for top 50% variable genes took a couple of hours
saveRDS(pkm, file="../int/wgcna/whi_preclustering.rds")


# Run WGCNA on full matrix of WHI beta values
system.time(bwnet <- blockwiseModules(t(betas_whi),
                                      blocks=pkm$clusters,
                                      power=6,
                                      networkType="unsigned",
                                      saveTOMs=T,
                                      saveTOMFileBase="../int/wgcna/whiTOM",
                                      verbose=2))
## blockwise network construction took 4.5 hrs when given pkm clustering blocks
saveRDS(bwnet, file="../int/wgcna/whi_net.rds")


# Calculate module eigengenes
system.time(MEs <- moduleEigengenes(t(betas_whi), bwnet$colors, excludeGrey=T))
save(MEs, file="../int/wgcna/whi_MEs.rds")
## eigengene calculation took 0.5 hrs


# Module preservation in FHS
multi_expr <- list(whi=list(data=t(betas_whi)),
                  fhs=list(data=t(betas_fhs[(
                    rownames(betas_fhs) %in% wgcnaInputCpgs), ])))
multi_color <- list(whi=whiNet$colors)
system.time(mp <- modulePreservation(multi_expr, 
                                     multi_color,
                                     networkType="unsigned",
                                     referenceNetworks=1,
                                     nPermutations=5,
                                     randomSeed=1,
                                     quickCor=1,
                                     verbose=3))
## modulePreservation takes ~40 mins. for full set of markers
saveRDS(mp, file="../int/wgcna/fhs_MP.rds")
