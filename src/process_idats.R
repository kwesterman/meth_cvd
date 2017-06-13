suppressMessages(silent <- lapply(c("tidyverse","minfi"), library, character.only=T))

args <- commandArgs(trailingOnly=T)
num_cores <- as.numeric(args[1])

## Sample sheet
targets <- read.metharray.sheet(base="../data/", pattern="^sample_sheet4.csv$") %>%
  dplyr::filter(Basename!="character(0)",  # Seems to appear in read.metharray.sheet output when a channel is missing
                Sample_Name!=21084)  # This sample is missing red channel .idat
save("targets", file="../int/targets.RData")

## Read in .idat intensity files
cl <- makePSOCKcluster(num_cores)
clusterExport(cl, "targets")  # Make targets data frame available to workers
silent <- clusterEvalQ(cl, library(minfi))  # Load minfi in all workers
print("Reading .idat files...")
rgSet_list <- parLapply(cl, splitIndices(nrow(targets), num_cores),  # Splits samples into equally-sized chunks
                        function(idx_set) read.metharray.exp(targets=targets[idx_set,]))
stopCluster(cl)
greens <- do.call(cbind, lapply(rgSet_list, getGreen))  # Extract greens from each subset
reds <- do.call(cbind, lapply(rgSet_list, getRed))  # Extract reds from each subset
rgSet <- RGChannelSet(Green=greens, Red=reds, annotation=annotation(rgSet_list[[1]]))  # Create full rgSet
rm(rgSet_list, greens, reds)  # Remove older objects to save memory
print("Saving rgSets...")
save("rgSet", file="../int/rgSet.RData")

## Generate p-values for probe detection
print("Running probe detection confidence tests...")
detP <- detectionP(rgSet)  # Generates a detection p-value for each probe in each sample
print("Saving detection p-values...")
save("detP", file="../int/detP.RData")


