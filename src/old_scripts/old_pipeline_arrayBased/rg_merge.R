suppressMessages(lapply(c("minfi"), library, character.only=T))

rgSet_files <- grep("rgSet[0-9]+.RData", list.files("../data/"), value=T)

rgSets <- lapply(rgSet_files, function(f) {
  print(f)
  load(paste0("../data/", f))
  rgSet
})
print(paste("All data subsets loaded. Memory used so far:", pryr::mem_used()))
stopifnot(all(unlist(lapply(rgSets, function(rgs) all(featureNames(rgs)==featureNames(rgSets[[1]]))))))  # Check that all CpG lengths and names are identical

## The merge step below requires on the order of 10MB per sample in available memory
grn <- do.call(cbind, lapply(rgSets, getGreen))
print("Green extracted.")
red <- do.call(cbind, lapply(rgSets, getRed))
print("Red extracted.")
print(paste("Memory used so far:", pryr::mem_used()))
rgSet <- RGChannelSet(Green=grn, Red=red)
load("../data/targets.RData")
annotation(rgSet) <- annotation(rgSets[[1]])
print(paste("rgSet created. Mem used so far:", pryr::mem_used()))

save("rgSet", file="../data/rgSet_full.RData")
