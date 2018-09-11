# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(
  c("tidyverse", "minfi", "RefFreeEWAS", "FlowSorted.Blood.450k", 
    "flashpcaR","sva"), library, character.only=T))

datasets <- commandArgs(TRUE)

estimate_cell_counts <- function(dataset, betas) {
  # Estimate cell counts using Houseman reference-based method
  print("...cell count estimation...")
  WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), 
                            rownames(betas)), ]
  stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
  est_cell_counts <- projectMix(WBC_subset, 
                                FlowSorted.Blood.450k.JaffeModelPars)
  est_cell_counts <- cbind(sampleKey=as.character(rownames(est_cell_counts)), 
                           data.frame(est_cell_counts))
  saveRDS(est_cell_counts, 
          file=paste0("../int/est_cell_counts_", dataset, ".rds"))
}

run_cpacor <- function(dataset) {
  # CPACOR (Lehne et al. 2015) to find control probe PCs that adjust for 
  # technical confounding
  print("...CPA method (control probe PCA)...")
  rg_set <- readRDS(paste0("../int/rg_set_", dataset, ".rds"))
  source("helpers.R")
  cpacor_fit <- run_CPA(rg_set)
  CP_PCs <- data.frame(cpacor_fit$x) %>%
    select(one_of(paste0("PC", 1:20))) %>%
    rename_all(funs(paste0("cp", .))) %>%  # To distinguish from standard PCA components
    rownames_to_column(var="sampleKey") %>%
    dplyr::filter(sampleKey %in% colnames(betas))  
      # B/c used rg_set as input, must trim samples that were filtered
  save("cpacor_fit", "CP_PCs", file=paste0("../int/CPACOR_", dataset, ".RData"))
}

run_pca <- function(dataset, betas, num_features=100000) {
  # PCA (using flashPCA method for speed)
  print(paste0("...PCA on ", num_features / 1000, "k most variable probes..."))
  probe_vars <- apply(betas, 1, var)
  high_var_indices <- order(probe_vars, decreasing=T)[1:num_features]
  high_var_mat <- betas[high_var_indices, ]
  pca_fit <- flashpca(t(high_var_mat), stand="sd", ndim=20)
  EVs <- setNames(data.frame(pca_fit$vectors), paste0("PC", 1:20))
  PCs <- cbind(sampleKey=as.character(colnames(high_var_mat)), EVs)
  save("pca_fit", "PCs", file=paste0("../int/pca_fit_", dataset, ".RData"))
}

run_combat <- function(dataset, betas) {
  if (is.null(combat_vars[[dataset]])) { 
    return(NULL)
  } else {
    print("...ComBat batch effect adjustment...")
    m_vals <- logit2(betas)
    metadata <- readRDS("../int/metadata.rds")
    metadata <- metadata[match(colnames(m_vals), metadata$sampleKey),]
    m_vals_combat <- ComBat(m_vals, batch=metadata[[combat_vars[[dataset]]]], 
                           mod=model.matrix(~event, data=metadata))
    betas_combat <- ilogit2(m_vals_combat)
    saveRDS(betas_combat, 
            file=paste0("../int/betas_qc_norm_filt_combat_", dataset, ".rds"), 
            compress=F)
  }
}

combat_vars <- list(fhs="center")

for (dataset in datasets) {
  print(paste("Dataset:", dataset))
  betas <- readRDS(paste0("../int/betas.qc.norm.filt_", dataset, ".rds"))
  silent <- estimate_cell_counts(dataset, betas)
  silent <- run_cpacor(dataset)
  silent <- run_pca(dataset, betas)
  silent <- run_combat(dataset, betas)
}

# if (!("betas.qc.norm.filt_combat.rds" %in% list.files("../int/"))) {
#   m_vals <- logit2(betas)
#   print("ComBat batch effect adjustment without outcome as covariate...")
#   m_vals_combat <- ComBat(m_vals, batch=metadata$study)
#   betas_combat <- ilogit2(m_vals_combat)
#   saveRDS(betas_combat, file="../int/betas.qc.norm.filt_combat.rds", compress=F)
#   print("PCA on basic ComBat-adjusted dataset...")
#   run_pca(betas_combat, 100000, "../int/PCA_combat.RData")
#   rm(m_vals_combat, betas_combat)
#   print("ComBat batch effect adjustment with outcome as covariate...")
#   m_vals_combat_events <- ComBat(m_vals, batch=metadata$study, mod=model.matrix(~event, data=metadata))
#   betas_combat_events <- ilogit2(m_vals_combat_events)
#   saveRDS(betas_combat_events, file="../int/betas.qc.norm.filt_combat_eventAdj.rds", compress=F)
#   print("PCA on ComBat-adjusted (with outcome as covariate) dataset...")
#   run_pca(betas_combat_events, 100000, "../int/PCA_combatWithOutcome.RData")
#   rm(m_vals, m_vals_combat_events, betas_combat_events)
