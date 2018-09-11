# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(c("tidyverse","minfi","gtools","sva","RefFreeEWAS","FlowSorted.Blood.450k","flashpcaR"), 
                                  library, character.only=T))

metaData <- readRDS("../int/metaData.rds")
betas <- readRDS("../int/betas.qc.norm.filt.rds")

## Ensure users for methylation and metadata match
metaData <- metaData[match(colnames(betas), metaData$sampleKey),]
stopifnot(all(metaData$sampleKey==colnames(betas)))

## Quick summary of M-value matrix
meth_matrix_dimensions <- dim(betas)
save("meth_matrix_dimensions", file="../int/meth_summary.RData")

## Estimate cell counts using Houseman reference-based method
if (!("estCellCounts.RData" %in% list.files("../int/"))) {
  print("Cell count estimation...")
  # betas <- ilogit2(Mvals)
  WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(betas)),]
  stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
  estCellCounts <- projectMix(WBC_subset, FlowSorted.Blood.450k.JaffeModelPars)
  estCellCounts <- cbind(sampleKey=as.character(rownames(estCellCounts)), data.frame(estCellCounts))
  save("estCellCounts", file="../int/estCellCounts.RData")  # Save as .RData with associated index no.
  print("Cell counts matrix saved.")
} else load("../int/estCellCounts.RData")


# ## SVA (Leek & Storey 2007) to find SVs capturing latent structure
# if (!("SVA.RData" %in% list.files("../int/"))) {
#   print("SVA...")
#   mod <- model.matrix(~event, data=metaData)
#   mod0 <- model.matrix(~1, data=metaData)
#   sva.fit <- sva(Mvals, mod, mod0, vfilter=5000)   # Raise this vfilter number?
#   SVs <- data.frame(sva.fit$sv) %>%
#     setNames(paste0("SV", 1:ncol(sva.fit$sv))) %>%
#     mutate(sampleKey=metaData$sampleKey)
#   save("sva.fit", "SVs", file="../int/SVA.RData")
#   print("Surrogate variables saved.")
# } else load("../int/SVA.RData")

## CPACOR (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
if (!("CPACOR.RData" %in% list.files("../int/"))) {
  print("CPA method (control probe PCA)...")
  rgSet <- readRDS("../int/rgSet.rds")
  source("helpers.R")
  cpacor.fit <- run_CPA(rgSet)
  CP_PCs <- data.frame(cpacor.fit$x) %>%
    select(one_of(paste0("PC",1:20))) %>%
    rename_all(paste0, "_cp") %>%  # To distinguish from standard PCA components
    rownames_to_column(var="sampleKey") %>%
    dplyr::filter(sampleKey %in% colnames(betas))  # B/c used rgSet as input, must trim samples that were filtered
  save("cpacor.fit", "CP_PCs", file="../int/CPACOR.RData")
  print("CPA results saved.")
  rm(rgSet)
} else load("../int/CPACOR.RData")

# ## ReFACTor (using their code with slight modifications)
# if (!("ReFACTor.RData" %in% list.files("../int/"))) {
#   print("ReFACTor...")
#   source("helpers.R")
#   rf_PCs <- refactor(ilogit2(Mvals), k=5)
#   colnames(rf_PCs) <- paste0(colnames(rf_PCs), "_rf")
#   rf_PCs <- cbind(sampleKey=colnames(Mvals), data.frame(rf_PCs))
#   save("rf_PCs", file="../int/ReFACTor.RData")
#   print("ReFACTor results saved.")
# } else load("../int/ReFACTor.RData")

run_pca <- function(dataMat, numFeatures, saveFileName) {
  probeVars <- apply(dataMat, 1, var)
  highVarIndices <- order(probeVars, decreasing=T)[1:numFeatures]
  highVarMat <- dataMat[highVarIndices,]
  pca.fit <- flashpca(t(highVarMat), stand="sd", ndim=20)
  EVs <- setNames(data.frame(pca.fit$vectors), paste0("PC", 1:20))
  PCs <- cbind(sampleKey=as.character(colnames(highVarMat)), EVs)
  save("pca.fit", "PCs", file=saveFileName)
}

## Basic PCA on most variable CpGs
if (!("PCA.RData" %in% list.files("../int/"))) {
  print("PCA on 100k most variable probes...")
  run_pca(betas, 100000, "../int/PCA.RData")
} else load("../int/PCA.RData")

if (!("betas.qc.norm.filt.combat.rds" %in% list.files("../int/"))) {
  Mvals <- logit2(betas)
  print("ComBat batch effect adjustment without outcome as covariate...")
  Mvals.combat <- ComBat(Mvals, batch=metaData$study)
  betas.combat <- ilogit2(Mvals.combat)
  saveRDS(betas.combat, file="../int/betas.qc.norm.filt.combat.rds", compress=F)
  print("PCA on basic ComBat-adjusted dataset...")
  run_pca(betas.combat, 100000, "../int/PCA_combat.RData")
  rm(Mvals.combat, betas.combat)
  print("ComBat batch effect adjustment with outcome as covariate...")
  Mvals.combat_events <- ComBat(Mvals, batch=metaData$study, mod=model.matrix(~event, data=metaData))
  betas.combat_events <- ilogit2(Mvals.combat_events)
  saveRDS(betas.combat_events, file="../int/betas.qc.norm.filt.combat_eventAdj.rds", compress=F)
  print("PCA on ComBat-adjusted (with outcome as covariate) dataset...")
  run_pca(betas.combat_events, 100000, "../int/PCA_combatWithOutcome.RData")
  rm(Mvals, Mvals.combat_events, betas.combat_events)
}



### Create master covariate/event data frame
nonMethData <- Reduce(function(x,y) inner_join(x, y, by="sampleKey"), 
                          list(metaData, estCellCounts, CP_PCs, PCs)) %>%
  distinct(subjID, .keep_all=T)  # Hacky...maybe revisit new sample sheet later or other way of dealing with duplicate shareids?
saveRDS(nonMethData, file="../int/nonMethData.rds")






