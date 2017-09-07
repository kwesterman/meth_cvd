suppressMessages(silent <- lapply(c("tidyverse","minfi","doParallel"), library, character.only=T))

## Sample sheet
methSheet <- read.metharray.sheet(base="../data/", pattern="^sample_sheet4.csv$")
phen <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv")
targets <- inner_join(methSheet, phen, by=c("Sample_Name"="shareid")) %>%
  dplyr::mutate(shareid=Sample_Name, 
                sex=c("M","F")[SEX], 
                age=as.numeric(AGE8)) %>%
  dplyr::select(shareid, Array, Slide, Basename, sex, age) %>%
  dplyr::filter(Basename!="character(0)",  # Seems to appear in read.metharray.sheet output when a channel is missing
                shareid!=21084)  # This sample is missing red channel .idat
save("targets", file="../int/targets.RData")

## Read in .idat intensity files
print("Reading .idat files...")
cl <- makePSOCKcluster(detectCores(), outfile="")
registerDoParallel(cl)
targetSplit <- splitIndices(nrow(targets), detectCores())  # Splits samples into equally-sized chunks
rgSet <- foreach(targetIdxSet=targetSplit, .packages="minfi",
                 .combine="combine", .multicombine=T) %dopar% 
  read.metharray.exp(targets=targets[targetIdxSet,])
stopCluster(cl)
# greens <- do.call(cbind, lapply(rgSet_list, getGreen))  # Extract greens from each subset
# reds <- do.call(cbind, lapply(rgSet_list, getRed))  # Extract reds from each subset
# rgSet <- RGChannelSet(Green=greens, Red=reds, annotation=annotation(rgSet_list[[1]]))  # Create full rgSet
# pData(rgSet) <- do.call(rbind, lapply(rgSet_list, pData))
# rm(rgSet_list, greens, reds)  # Remove older objects to save memory
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


