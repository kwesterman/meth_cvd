suppressMessages(silent <- lapply(c("tidyverse","minfi","doParallel"), library, character.only=T))

## Sample sheet
load("../int/sampleSheet.RData")
targets <- sampleSheet %>%
  dplyr::mutate(Basename=paste0("../data/raw_methylation/", sampleKey)) %>%
  dplyr::filter(file.exists(paste0(Basename, "_Red.idat")),
                file.exists(paste0(Basename, "_Grn.idat"))) %>%
  distinct()

## Read in .idat intensity files
print("Reading .idat files...")
# rgSet <- read.metharray.exp(targets=targets, verbose=T)
numCores <- min(16, detectCores())
cl <- makePSOCKcluster(numCores)
registerDoParallel(cl)
targetSplit <- splitIndices(nrow(targets), 2*numCores)  # Splits samples into equally-sized chunks
rgSet <- foreach(targetIdxSet=targetSplit, .packages="minfi", 
                      .combine=BiocGenerics::combine, .multicombine=T) %dopar%
  read.metharray.exp(targets=targets[targetIdxSet,])
stopCluster(cl)

print("Saving rgSet...")
save("rgSet", file="../int/rgSet.RData")

## Generate p-values for probe detection
print("Running probe detection confidence tests...")
detP <- detectionP(rgSet)  # Generates a detection p-value for each probe in each sample
print("Saving detection p-values...")
save("detP", file="../int/detP.RData")

## Produce U vs. M intensity plot for visual QC
print("Creating QC intensity plot...")
jpeg("../output/qcIntensityPlot.jpg")  # Plot of Unmeth vs. Meth median intensities
plotQC(getQC(preprocessRaw(rgSet)))
dev.off()


