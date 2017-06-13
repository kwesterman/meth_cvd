
lapply(c("knitr","tidyverse"), library, character.only=T)

load("../data/methylation_covariates.Rdata")  # Technical variables (from sample sheet) and cell counts
phen2 <- read.csv("../data/Off_Exam8_ex1_phen2.csv")  # Has BMI
phenAll <- read.csv("../data/Off_Exam8_phen_cov_all.csv")  # Has age and smoking status

covData <- phenAll %>%
  dplyr::select(shareid, SEX, AGE8, smoking_now) %>%
  inner_join(phen2, by="shareid") %>%
  dplyr::select(shareid, dbGaP_Subject_ID, SEX, AGE8, smoking_now, bmi) %>%
  dplyr::rename(sex=SEX, age=AGE8) %>%
  inner_join(meth_covData, by=c("shareid"="Sample_Name"))

save("covData", file="../data/covariates.Rdata")


