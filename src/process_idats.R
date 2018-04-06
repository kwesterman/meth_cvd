suppressMessages(silent <- lapply(c("tidyverse","minfi","doParallel"), library, character.only=T))

datasets <- commandArgs(TRUE)  # Create list of datasets to use from command-line arguments

create_targets <- function(dataset) {
  # For a given dataset, use the sample sheet to generate a "targets" data frame to use in .idat read-in
  sampleSheet <- readRDS("../int/sampleSheet.rds")
  targets <- sampleSheet %>%
    dplyr::mutate(Basename=paste0("../data/", dataset, "/meth/", sampleKey)) %>%
    dplyr::filter(study==dataset) %>%
    dplyr::filter(file.exists(paste0(Basename, "_Red.idat")),
                  file.exists(paste0(Basename, "_Grn.idat"))) %>%
    distinct()
  targets
}

read_idats <- function(dataset, targets, numCores) {
  print(paste0("...reading .idat files..."))
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  targetSplit <- splitIndices(nrow(targets), 2*numCores)  # Splits samples into equally-sized chunks
  rgSet <- foreach(targetIdxSet=targetSplit, .packages="minfi", 
                   .combine=BiocGenerics::combine, .multicombine=T) %dopar%
    read.metharray.exp(targets=targets[targetIdxSet,])
  stopCluster(cl)
  
  print("...saving rgSet...")
  saveRDS(rgSet, file=paste0("../int/rgSet_", dataset, ".rds"), compress=F)
  
  rgSet
}

run_detectionP <- function(dataset, rgSet) {
  # Generate p-values for probe detection
  print("...running probe detection confidence tests...")
  detP <- detectionP(rgSet)  # Generates a detection p-value for each probe in each sample
  print("...saving detection p-values...")
  saveRDS(detP, file=paste0("../int/detP_", dataset, ".rds"), compress=F)
}

create_qcIntensityPlot <- function(dataset, rgSet) {
  # Produce U vs. M intensity plot for visual QC
  print("...creating QC intensity plot...")
  jpeg(paste0("../output/qcIntensityPlot_", dataset, ".jpg"))  # Plot of Unmeth vs. Meth median intensities
  plotQC(getQC(preprocessRaw(rgSet)))
  dev.off()
}

for (dataset in datasets) {
  print(paste0("Dataset: ", dataset))
  targetsDF <- create_targets(dataset)
  rgSet <- read_idats(dataset, targetsDF, 16)
  silent <- run_detectionP(dataset, rgSet)
  silent <- create_qcIntensityPlot(dataset, rgSet)
}


