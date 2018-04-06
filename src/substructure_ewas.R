# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(c("tidyverse","minfi","RefFreeEWAS","FlowSorted.Blood.450k","flashpcaR","sva"), 
                                  library, character.only=T))

datasets <- commandArgs(TRUE)

estimate_cell_counts <- function(dataset, betas) {
  # Estimate cell counts using Houseman reference-based method
  print("...cell count estimation...")
  WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(betas)),]
  stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
  estCellCounts <- projectMix(WBC_subset, FlowSorted.Blood.450k.JaffeModelPars)
  estCellCounts <- cbind(sampleKey=as.character(rownames(estCellCounts)), data.frame(estCellCounts))
  saveRDS(estCellCounts, file=paste0("../int/estCellCounts_", dataset, ".rds"))
}

run_cpacor <- function(dataset) {
  # CPACOR (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
  print("...CPA method (control probe PCA)...")
  rgSet <- readRDS(paste0("../int/rgSet_", dataset, ".rds"))
  source("helpers.R")
  cpacor.fit <- run_CPA(rgSet)
  CP_PCs <- data.frame(cpacor.fit$x) %>%
    select(one_of(paste0("PC",1:20))) %>%
    rename_all(funs(paste0("cp",.))) %>%  # To distinguish from standard PCA components
    rownames_to_column(var="sampleKey") %>%
    dplyr::filter(sampleKey %in% colnames(betas))  # B/c used rgSet as input, must trim samples that were filtered
  save("cpacor.fit", "CP_PCs", file=paste0("../int/CPACOR_", dataset, ".RData"))
}

run_pca <- function(dataset, betas, numFeatures=100000) {
  # PCA (using flashPCA method for speed)
  print(paste0("...PCA on ", numFeatures/1000, "k most variable probes..."))
  probeVars <- apply(betas, 1, var)
  highVarIndices <- order(probeVars, decreasing=T)[1:numFeatures]
  highVarMat <- betas[highVarIndices,]
  pca.fit <- flashpca(t(highVarMat), stand="sd", ndim=20)
  EVs <- setNames(data.frame(pca.fit$vectors), paste0("PC", 1:20))
  PCs <- cbind(sampleKey=as.character(colnames(highVarMat)), EVs)
  save("pca.fit", "PCs", file=paste0("../int/PCA.fit_", dataset, ".RData"))
}

run_combat <- function(dataset, betas) {
  if (is.null(combat_vars[[dataset]])) return(NULL)
  else {
    print("...ComBat batch effect adjustment...")
    Mvals <- logit2(betas)
    metaData <- readRDS("../int/metaData.rds")
    metaData <- metaData[match(colnames(Mvals), metaData$sampleKey),]
    Mvals.combat <- ComBat(Mvals, batch=metaData[[combat_vars[[dataset]]]], 
                           mod=model.matrix(~event, data=metaData))
    betas.combat <- ilogit2(Mvals.combat)
    saveRDS(betas.combat, file=paste0("../int/betas.qc.norm.filt.combat_", dataset, ".rds"), compress=F)
  }
}

combat_vars <- list(whi=NULL,
                    fhs="center")

for (dataset in datasets) {
  print(paste("Dataset:", dataset))
  betas <- readRDS(paste0("../int/betas.qc.norm.filt_", dataset, ".rds"))
  silent <- estimate_cell_counts(dataset, betas)
  silent <- run_cpacor(dataset)
  silent <- run_pca(dataset, betas)
  silent <- run_combat(dataset, betas)
}

# if (!("betas.qc.norm.filt.combat.rds" %in% list.files("../int/"))) {
#   Mvals <- logit2(betas)
#   print("ComBat batch effect adjustment without outcome as covariate...")
#   Mvals.combat <- ComBat(Mvals, batch=metaData$study)
#   betas.combat <- ilogit2(Mvals.combat)
#   saveRDS(betas.combat, file="../int/betas.qc.norm.filt.combat.rds", compress=F)
#   print("PCA on basic ComBat-adjusted dataset...")
#   run_pca(betas.combat, 100000, "../int/PCA_combat.RData")
#   rm(Mvals.combat, betas.combat)
#   print("ComBat batch effect adjustment with outcome as covariate...")
#   Mvals.combat_events <- ComBat(Mvals, batch=metaData$study, mod=model.matrix(~event, data=metaData))
#   betas.combat_events <- ilogit2(Mvals.combat_events)
#   saveRDS(betas.combat_events, file="../int/betas.qc.norm.filt.combat_eventAdj.rds", compress=F)
#   print("PCA on ComBat-adjusted (with outcome as covariate) dataset...")
#   run_pca(betas.combat_events, 100000, "../int/PCA_combatWithOutcome.RData")
#   rm(Mvals, Mvals.combat_events, betas.combat_events)
