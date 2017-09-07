suppressMessages(lapply(c("dplyr","minfi"), library, character.only=T))

dataDir <- "../data/"  # At the moment, this only contains a few samples for development
targets <- read.metharray.sheet(base=dataDir, pattern="^sample_sheet4.csv$") %>%
  dplyr::filter(Basename!="character(0)",
                Sample_Name!=21084)

save("targets", file="../data/targets.RData")
