library(tidyverse)
library(WGCNA)
library(survival)

# betas_topVar <- betas[beta_variances>quantile(beta_variances,0.5),]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), stringsAsFactors=F)

tssNoIslands <- anno450k %>% 
  dplyr::rename(cpg=Name, group=UCSC_RefGene_Group) %>%
  dplyr::select(cpg, group, Relation_to_Island) %>%
  separate_rows(group, sep=";") %>% 
  filter(group %in% c("TSS200","TSS1500"), 
         Relation_to_Island!="Island")


enableWGCNAThreads(2)
# use enableWGCNAThreads() instead of allow...()
# why/when does pickSoftThreshold func slow down? matter how many samples?

# powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# system.time(sft <- pickSoftThreshold(methDat_subset, powerVector = powers, verbose = 5))
## soft threshold choice on 500 samples for top 50% variable genes took 5.8 hrs --> indicates that power ~ 8 is good

# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## based on visual analysis, preliminary choice of beta=8

betas <- readRDS("../int/betas.qc.norm.filt_whi.rds")
# beta_variances <- apply(betas, 1, var)
# betas_topVar <- betas[beta_variances>quantile(beta_variances,0.75),]
methDat <- t(betas[rownames(betas) %in% tssNoIslands$cpg,])
bwnet <- blockwiseModules(methDat,
                          maxBlockSize=20000,
                          power=8,
                          networkType="signed",
                          deepSplit=2,
                          verbose=2)
saveRDS(bwnet, "../int/wgcna_network_whi.rds")

MEs0 <- moduleEigengenes(methDat, bwnet$colors)$eigengenes
MEs <- orderMEs(MEs0)
saveRDS(MEs, "../int/wgcna_eigengenes_whi.rds")



multiExpr <- 




whiExprSmall <- t(betas_whi[rownames(betas_whi) %in% tssNoIslands$cpg,])
fhsExprSmall <- t(betas_fhs[rownames(betas_fhs) %in% tssNoIslands$cpg,])
whiColors <- net$colors
multiExpr <- list(whi=list(data=whiExprSmall), fhs=list(data=fhsExprSmall))
multiColor <- list(whi=whiColors)
system.time(mp <- modulePreservation(multiExpr, multiColor,
                                     networkType="unsigned",
                                     referenceNetworks=1,
                                     nPermutations=50,
                                     randomSeed=1,
                                     quickCor=1,
                                     verbose=3))




# betas <- readRDS("../int/betas.qc.norm.filt_fhs.rds")
# beta_variances <- apply(betas, 1, var)
# betas_topVar <- betas[beta_variances>quantile(beta_variances,0.75),]
# methDat <- t(betas_topVar)
# bwnet <- blockwiseModules(methDat,
#                           maxBlockSize=20000,
#                           power=8,
#                           networkType="unsigned",
#                           deepSplit=2,
#                           verbose=2)
# saveRDS(bwnet, "../int/wgcna_network_fhs.rds")
# 
# MEs0 <- moduleEigengenes(methDat, bwnet$colors)$eigengenes
# MEs <- orderMEs(MEs0)
# saveRDS(MEs, "../int/wgcna_eigengenes_fhs.rds")