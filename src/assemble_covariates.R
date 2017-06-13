# Prepare covariate data for EWAS analysis of incident CVD events

suppressMessages(silent <- lapply(c("tidyverse","minfi","sva","RefFreeEWAS","FlowSorted.Blood.450k"), 
                                  library, character.only=T))

### PHENOTYPIC COVARIATES
print("Phenotype variables...")
# Load metadata/phenotype data   ##### MASKING OF DPLYR ISSUE
phenAll <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
  dplyr::rename(sex=SEX, age=AGE8) %>%
  dplyr::select(shareid, dbGaP_Subject_ID, sex, age, smoking_now, CVD_med, HT_med, T2D_med, SBP)
phen2 <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
  dplyr::select(shareid, bmi)
labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
  dplyr::select(-dbGaP_Subject_ID, -IDTYPE)

phenoData <- phenAll %>%
  inner_join(phen2, by="shareid") %>%
  inner_join(labValData, by="shareid")
save("phenoData", file="../int/phenoData.RData")

### CVD DATA
print("Cardiovascular event data...")
# Load and clean event data
soe <- read.delim("../data/events/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)

eventData <- soe %>%
  inner_join(labValData, by="shareid") %>%
  dplyr::mutate(timeToEvent=DATE-drawdate) %>%  # Days between sample and event (NOTE: specific values here will change)
  dplyr::filter(timeToEvent>0) %>%  # Don't keep events that occurred prior to methylation collection
  dplyr::filter(EVENT %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
  group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
  dplyr::slice(which.min(timeToEvent)) %>%  # Only keep first relevant event per person -- allow people with existing CVD for now
  ungroup() %>%
  dplyr::select(shareid, EVENT, drawdate, timeToEvent)
save("eventData", file="../int/eventData.RData")


### MOLECULAR COVARIATES
print("Molecular covariates:")
load("../int/Mvals.RData")
betas <- ilogit2(Mvals)
sampleData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
  dplyr::mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_")) %>%
  dplyr::rename(shareid=Sample_Name) %>%
  dplyr::select(shareid, sampleKey, Sentrix_ID, Sentrix_Position) %>%
  dplyr::slice(match(colnames(Mvals), sampleKey))
save("sampleData", file="../int/sampleData.RData")

## Estimate cell counts using Houseman reference-based method
print("Cell count estimation...")
WBC_subset <- betas[match(rownames(FlowSorted.Blood.450k.JaffeModelPars), rownames(betas)),]
stopifnot(nrow(WBC_subset)==nrow(FlowSorted.Blood.450k.JaffeModelPars))
estCellCounts <- projectMix(WBC_subset, FlowSorted.Blood.450k.JaffeModelPars)
save("estCellCounts", file="../int/estCellCounts.RData")  # Save as .RData with associated index no.
print("Cell counts matrix saved.")

## SVA (Leek & Storey 2007) to find SVs capturing latent structure
print("SVA...")
svaData <- mutate(sampleData, incidentCVD=ifelse(shareid %in% eventData$shareid, 1, 0))
mod <- model.matrix(~as.factor(incidentCVD), data=svaData)
mod0 <- model.matrix(~1, data=svaData)
cleanMvals <- Mvals[apply(Mvals, 1, function(r) all(!is.na(r) & is.finite(r))),]
svobj <- sva(cleanMvals, mod, mod0, vfilter=5000)   # Raise this vfilter number?
SVs <- svobj$sv
save("SVs", file="../int/SVs.RData")
print("Surrogate variables saved.")
rm(Mvals, betas)

## CPA (Lehne et al. 2015) to find control probe PCs that adjust for technical confounding
print("CPA method (control probe PCA)...")
load("../int/rgSet.RData")
source("Lehne2015.R")
CP_PCs <- run_CPA(rgSet)
save("CP_PCs", file="../int/CP_PCs.RData")
print("CPA results saved.")


# ######### IF I WANT TO MERGE INTO A MASTER COVARIATE/EVENT MATRIX ##########
# # Load technical sample covariates and merge with other sample data
# technicalData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
#   mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_")) %>%
#   rename(shareid=Sample_Name) %>%
#   select(shareid, sampleKey, Sentrix_ID, Sentrix_Position)
# 
# sampleData <- technicalData %>%
#   inner_join(estCellCounts, by="sampleKey") %>%  # Add estimated cell counts
#   inner_join(CP_PCs, by="sampleKey")  # Add CPA control probe PCs
# save("sampleData", file="../int/sampleData.RData")





