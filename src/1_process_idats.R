suppressMessages(silent <- lapply(c("tidyverse", "minfi", "doParallel"), library, character.only=T))

datasets <- commandArgs(TRUE)  # Create list of datasets to use from command-line arguments

create_targets <- function(dataset) {
  # For a given dataset, use the sample sheet to generate a 
  # "targets" data frame to use in .idat read-in
  sample_sheet <- readRDS("../int/sample_sheet.rds")
  targets <- sample_sheet %>%
    dplyr::mutate(Basename=paste0("../data/", dataset, 
                                  "/meth/", sampleKey)) %>%
    dplyr::filter(study == dataset) %>%
    dplyr::filter(file.exists(paste0(Basename, "_Red.idat")),
                  file.exists(paste0(Basename, "_Grn.idat"))) %>%
    distinct()
  targets
}

read_idats <- function(dataset, targets, numCores) {
  print(paste0("...reading .idat files..."))
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  target_split <- splitIndices(nrow(targets), 2 * numCores)  # Splits samples into equally-sized chunks
  rg_set <- foreach(target_idx_set=target_split, .packages="minfi", 
                   .combine=BiocGenerics::combine, .multicombine=T) %dopar%
    read.metharray.exp(targets=targets[target_idx_set, ])
  stopCluster(cl)
  
  print("...saving rg_set...")
  saveRDS(rg_set, file=paste0("../int/rg_set_", dataset, ".rds"), compress=F)
  
  rg_set
}

run_detectionP <- function(dataset, rg_set) {
  # Generate p-values for probe detection
  print("...running probe detection confidence tests...")
  detP <- detectionP(rg_set)  # Generates a detection p-value for each probe in each sample
  print("...saving detection p-values...")
  saveRDS(detP, file=paste0("../int/detP_", dataset, ".rds"), compress=F)
}

create_qc_intensity_plot <- function(dataset, rg_set) {
  # Produce U vs. M intensity plot for visual QC
  print("...creating QC intensity plot...")
  jpeg(paste0("../output/qcIntensityPlot_", dataset, ".jpg"))  # Plot of Unmeth vs. Meth median intensities
  plotQC(getQC(preprocessRaw(rg_set)))
  dev.off()
}

for (dataset in datasets) {
  print(paste0("Dataset: ", dataset))
  targets <- create_targets(dataset)
  rg_set <- read_idats(dataset, targets, 16)
  silent <- run_detectionP(dataset, rg_set)
  silent <- create_qc_intensity_plot(dataset, rg_set)
}


