suppressMessages(silent <- lapply(c("tidyverse","minfi","doParallel"), library, character.only=T))

## Sample sheet
sampleSheet <- readRDS("../int/sampleSheet.rds")
targets <- sampleSheet %>%
  dplyr::mutate(Basename=paste0("../data/", study, "/meth/", sampleKey)) %>%
  dplyr::filter(file.exists(paste0(Basename, "_Red.idat")),
                file.exists(paste0(Basename, "_Grn.idat"))) %>%
  distinct()

## Read in .idat intensity files
print("Reading .idat files...")

numCores <- 16
cl <- makePSOCKcluster(numCores)
registerDoParallel(cl)
targetSplit <- splitIndices(nrow(targets), 2*numCores)  # Splits samples into equally-sized chunks
rgSet <- foreach(targetIdxSet=targetSplit, .packages="minfi", 
                      .combine=BiocGenerics::combine, .multicombine=T) %dopar%
  read.metharray.exp(targets=targets[targetIdxSet,])
stopCluster(cl)

print("Saving rgSet...")
saveRDS(rgSet, file="../int/rgSet.rds")

## Generate p-values for probe detection
print("Running probe detection confidence tests...")
detP <- detectionP(rgSet)  # Generates a detection p-value for each probe in each sample
print("Saving detection p-values...")
saveRDS(detP, file="../int/detP.rds", compress=F)

## Produce U vs. M intensity plot for visual QC
print("Creating QC intensity plot...")
jpeg("../output/qcIntensityPlot.jpg")  # Plot of Unmeth vs. Meth median intensities
plotQC(getQC(preprocessRaw(rgSet)))
dev.off()


