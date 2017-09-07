# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(c("tidyverse","minfi","gtools","sva","RefFreeEWAS","FlowSorted.Blood.450k"), 
                                  library, character.only=T))

load("../int/metaData.RData")
load("../int/Mvals.RData")

## Ensure users for methylation and metadata match
metaData <- metaData[match(colnames(Mvals), metaData$sampleKey),]
stopifnot(all(metaData$sampleKey==colnames(Mvals)))

## Quick summary of M-value matrix
meth_matrix_dimensions <- dim(Mvals)
save("meth_matrix_dimensions", file="meth_summary.RData")

## Estimate cell counts using Houseman reference-based method
if (!("estCellCounts.RData" %in% list.files("../int/"))) {
  print("Cell count estimation...")
  betas <- ilogit2(Mvals)
  WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(betas)),]
  stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
  estCellCounts <- projectMix(WBC_subset, FlowSorted.Blood.450k.JaffeModelPars)
  estCellCounts <- cbind(sampleKey=as.character(rownames(estCellCounts)), data.frame(estCellCounts))
  save("estCellCounts", file="../int/estCellCounts.RData")  # Save as .RData with associated index no.
  print("Cell counts matrix saved.")
  rm(betas)
} else load("../int/estCellCounts.RData")


## SVA (Leek & Storey 2007) to find SVs capturing latent structure
if (!("SVA.RData" %in% list.files("../int/"))) {
  print("SVA...")
  mod <- model.matrix(~event, data=metaData)
  mod0 <- model.matrix(~1, data=metaData)
  sva.fit <- sva(Mvals, mod, mod0, vfilter=5000)   # Raise this vfilter number?
  SVs <- data.frame(sva.fit$sv) %>%
    setNames(paste0("SV", 1:ncol(sva.fit$sv))) %>%
    mutate(sampleKey=metaData$sampleKey)
  save("sva.fit", "SVs", file="../int/SVA.RData")
  print("Surrogate variables saved.")
} else load("../int/SVA.RData")

## CPACOR (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
if (!("CPACOR.RData" %in% list.files("../int/"))) {
  print("CPA method (control probe PCA)...")
  load("../int/rgSet.RData")
  source("helpers.R")
  cpacor.fit <- run_CPA(rgSet)
  CP_PCs <- cpacor.fit$x
  CP_PCs <- data.frame(cpacor.fit$x) %>%
    rename_all(paste0, "_cp") %>%  # To distinguish from standard PCA components
    rownames_to_column(var="sampleKey") %>%
    dplyr::filter(sampleKey %in% colnames(Mvals))  # B/c used rgSet as input, must trim samples that were filtered
  save("cpacor.fit", "CP_PCs", file="../int/CPACOR.RData")
  print("CPA results saved.")
  rm(rgSet)
} else load("../int/CPACOR.RData")

## ReFACTor (using their code with slight modifications)
if (!("ReFACTor.RData" %in% list.files("../int/"))) {
  print("ReFACTor...")
  source("helpers.R")
  rf_PCs <- refactor(ilogit2(Mvals), k=5)
  colnames(rf_PCs) <- paste0(colnames(rf_PCs), "_rf")
  rf_PCs <- cbind(sampleKey=colnames(Mvals), data.frame(rf_PCs))
  save("rf_PCs", file="../int/ReFACTor.RData")
  print("ReFACTor results saved.")
} else load("../int/ReFACTor.RData")

## Basic PCA on most variable CpGs
if (!("PCA.RData" %in% list.files("../int/"))) {
  print("PCA on 50k most variable probes...")
  probeVars <- apply(Mvals, 1, var) ### CLEANING NECESSARY???
  highVarIndices <- order(probeVars, decreasing=T)[1:50000]
  highVarProbes <- Mvals[highVarIndices,]
  pca.fit <- prcomp(t(highVarProbes))
  jpeg("../output/pca_screePlot.jpg")
  screeplot(pca.fit, npcs=25, type="lines")
  dev.off()
  PCs <- pca.fit$x[,1:25]
  PCs <- cbind(sampleKey=as.character(colnames(highVarProbes)), data.frame(PCs))
  save("pca.fit", "PCs", file="../int/PCA.RData")
} else load("../int/PCA.RData")

### Create master covariate/event data frame
nonMethData <- Reduce(function(x,y) inner_join(x, y, by="sampleKey"), 
                          list(metaData, estCellCounts, SVs, CP_PCs, rf_PCs, PCs)) %>%
  distinct(shareid, .keep_all=T)  # Hacky...maybe revisit new sample sheet later or other way of dealing with duplicate shareids?
save("nonMethData", file="../int/nonMethData.RData")






