suppressMessages(silent <- lapply(c("tidyverse"), library, character.only=T))

# Load metadata/phenotype data
phenAll <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
  rename(sex=SEX, age=AGE8) %>%
  select(shareid, dbGaP_Subject_ID, sex, age, smoking_now, CVD_med, HT_med, T2D_med, SBP)
phen2 <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
  select(shareid, bmi)
labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
  select(-dbGaP_Subject_ID, -IDTYPE)

phenoData <- phenAll %>%
  inner_join(phen2, by="shareid") %>%
  inner_join(labValData, by="shareid")
save("phenoData", file="../int/phenoData.RData")

# Load control probe PCs
load("../int/CP_PCs.RData")
CP_PCs <- data.frame(CP_PCs) %>%
  rownames_to_column(var="sampleKey") %>%
  select(sampleKey, one_of(paste0("PC",1:30,"_cp")))

# Load cell counts
estCellCount_files <- grep("cellCounts_[0-9]+.RData", list.files("../int/"), value=T)
estCellCount_list <- lapply(paste0("../int/", estCellCount_files), function(x) {load(x); cellCounts})
estCellCounts <- do.call(rbind, estCellCount_list) %>%
  data.frame() %>%
  rownames_to_column(var="sampleKey")

# Load technical sample covariates and merge with other sample data
technicalData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
  mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_")) %>%
  rename(shareid=Sample_Name) %>%
  select(shareid, sampleKey, Sentrix_ID, Sentrix_Position)

sampleData <- technicalData %>%
  inner_join(estCellCounts, by="sampleKey") %>%  # Add estimated cell counts
  inner_join(CP_PCs, by="sampleKey")  # Add CPA control probe PCs
save("sampleData", file="../int/sampleData.RData")


