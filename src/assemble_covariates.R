# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(c("tidyverse","minfi","gtools","sva","RefFreeEWAS","FlowSorted.Blood.450k"), 
                                  library, character.only=T))

### PHENOTYPIC COVARIATES
print("Phenotype variables...")
# Load metadata/phenotype data
phenAll <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
  dplyr::rename(sex=SEX, age=AGE8) %>%
  dplyr::select(shareid, dbGaP_Subject_ID, sex, age, smoking_now, CVD_med, HT_med, T2D_med, SBP)
phen2 <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
  dplyr::select(shareid, bmi)
labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
  dplyr::select(-dbGaP_Subject_ID, -IDTYPE)

phenoData <- phenAll %>%
  left_join(phen2, by="shareid") %>%
  left_join(labValData, by="shareid")
save("phenoData", file="../int/phenoData.RData")


### SAMPLE DATA
load("../int/Mvals.RData")
cleanMvals <- Mvals[apply(Mvals, 1, function(r) all(!is.na(r) & is.finite(r))),]

print("Sample metadata...")
sampleData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
  dplyr::mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_"),
                Sentrix_Row=substr(Sentrix_Position, 1, 3),
                Sentrix_Col=substr(Sentrix_Position, 4, 6)) %>%
  dplyr::rename(shareid=Sample_Name) %>%
  dplyr::select(shareid, sampleKey, Sentrix_ID, Sentrix_Row, Sentrix_Col) %>%
  dplyr::slice(match(colnames(cleanMvals), sampleKey))
save("sampleData", file="../int/sampleData.RData")

### CVD DATA
print("Cardiovascular event data...")
# Load and clean event data
soe <- read.delim("../data/events/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)

eventData <- soe %>%
  inner_join(labValData, by="shareid") %>%
  dplyr::mutate(timeToEvent=DATE-drawdate) %>%  # Days between sample and event (NOTE: specific values here will change)
  dplyr::filter(EVENT %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
  dplyr::filter(timeToEvent>0) %>%  # Don't keep events that occurred prior to methylation collection
  group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
  dplyr::slice(which.min(timeToEvent)) %>%  # Only keep first relevant event per person -- allow people with existing CVD for now
  ungroup() %>%
  dplyr::select(shareid, EVENT, drawdate, timeToEvent)
save("eventData", file="../int/eventData.RData")  # This is just the events, so not all subjects included

survData <- sampleData %>%  # Basis for including all relevant samples (not just those with events)
  left_join(eventData, by="shareid") %>% 
  left_join(labValData, by="shareid") %>%
  dplyr::mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
                timeToEvent=na.replace(timeToEvent, 1000)) %>%  ## NOTE: this 800 needs to be switched for an *actual* number
  dplyr::select(shareid, event, timeToEvent)
  
### TECHNICAL COVARIATES

## Estimate cell counts using Houseman reference-based method
if (!("estCellCounts.RData" %in% list.files("../int/"))) {
  print("Cell count estimation...")
  betas <- ilogit2(Mvals)
  WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(betas)),]
  stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
  estCellCounts <- projectMix(WBC_subset, FlowSorted.Blood.450k.JaffeModelPars)
  estCellCounts <- cbind(sampleKey=as.character(rownames(estCellCounts)), data.frame(estCellCounts))
  save("estCellCounts", file="../int/estCellCounts.RData")  # Save as .RData with associated index no.
  print("Cell counts matrix saved.")
  rm(betas)
} else load("../int/estCellCounts.RData")


## SVA (Leek & Storey 2007) to find SVs capturing latent structure
if (!("SVs.RData" %in% list.files("../int/"))) {
  print("SVA...")
  svaData <- mutate(sampleData, incidentCVD=ifelse(shareid %in% eventData$shareid, 1, 0))
  mod <- model.matrix(~as.factor(incidentCVD), data=svaData)
  mod0 <- model.matrix(~1, data=svaData)
  svobj <- sva(cleanMvals, mod, mod0, vfilter=5000)   # Raise this vfilter number?
  SVs <- svobj$sv
  colnames(SVs) <- paste0("SV", 1:ncol(SVs))
  SVs <- cbind(sampleKey=as.character(svaData$sampleKey), data.frame(SVs))
  save("SVs", file="../int/SVs.RData")
  print("Surrogate variables saved.")
} else load("../int/SVs.RData")

## CPA (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
if (!("CP_PCs.RData" %in% list.files("../int/"))) {
  print("CPA method (control probe PCA)...")
  load("../int/rgSet.RData")
  source("helpers.R")
  CP_PCs <- run_CPA(rgSet)
  CP_PCs <- CP_PCs[rownames(estCellCounts),]  # B/c used rgSet as input, must trim samples that were filtered
  CP_PCs <- cbind(sampleKey=as.character(rownames(CP_PCs)), data.frame(CP_PCs))
  save("CP_PCs", file="../int/CP_PCs.RData")
  print("CPA results saved.")
  rm(rgSet)
} else load("../int/CP_PCs.RData")

## ReFACTor (using their code with slight modifications)
if (!("rf_PCs.RData" %in% list.files("../int/"))) {
  print("ReFACTor...")
  source("helpers.R")
  rf_PCs <- refactor(ilogit2(Mvals), k=5)
  names(rf_PCs) <- paste(names(rf_PCs), "rf", sep="_")
  rf_PCs <- cbind(sampleKey=colnames(Mvals), data.frame(rf_PCs))
  save("rf_PCs", file="../int/rf_PCs.RData")
  print("ReFACTor results saved.")
} else load("../int/rf_PCs.RData")

## Basic PCA on most variable CpGs
if (!("PCs.RData" %in% list.files("../int/"))) {
  print("PCA on 50k most variable probes...")
  probeVars <- apply(cleanMvals, 1, var)
  highVarIndices <- order(probeVars, decreasing=T)[1:50000]
  highVarProbes <- cleanMvals[highVarIndices,]
  
  prcomp.fit <- prcomp(t(highVarProbes))
  PCs <- prcomp.fit$x[,1:20]
  PCs <- cbind(sampleKey=as.character(colnames(highVarProbes)), data.frame(PCs))
  save("PCs", file="../int/PCs.RData")
} else load("../int/PCs.RData")

## Methylation-based sample similarity matrix (most variable CpGs)
if (!("Kmat.RData" %in% list.files("../int/"))) {
  probeVars <- apply(cleanMvals, 1, var)
  highVarIndices <- order(probeVars, decreasing=T)[1:20000]
  highVarProbes <- cleanMvals[highVarIndices,]
  Kmat <- t(highVarProbes) %*% highVarProbes / nrow(highVarProbes)
  save("Kmat", file="../int/Kmat.RData")
}

technicalCovars <- Reduce(function(x,y) inner_join(x, y, by="sampleKey"), 
                          list(sampleData, estCellCounts, SVs, CP_PCs, rf_PCs, PCs))

### Create master covariate/event data frame

nonMethData <- technicalCovars %>%
  distinct(shareid, .keep_all=T) %>%  # Hacky...maybe revisit new sample sheet later or other way of dealing with duplicate shareids?
  inner_join(phenoData, by="shareid") %>%
  inner_join(survData, by="shareid")
save("nonMethData", file="../int/nonMethData.RData")

## Investigate correlations of covariates with top 10 PCs

PC_df <- inner_join(PCs[,1:11], 
                    nonMethData[,unique(c("sampleKey","sex","age","event",names(estCellCounts)))],
                    by="sampleKey") %>%
  dplyr::select(-sampleKey)
PC_corMat <- cor(PC_df)[-(1:10),]
save("PC_corMat", file="../output/PC_corMat.RData")
jpeg("../output/PC_heatmap.jpg"); heatmap(PC_corMat, Rowv=NA); dev.off()




