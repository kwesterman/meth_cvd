# Assemble metadata (demographics, CVD events, lab values, etc.) on subjects

suppressMessages(silent <- lapply(c("tidyverse"), library, character.only=T))

### SAMPLE DATA
print("Sample metadata...")

sampleData <- read.csv("../data/raw_methylation/sample_sheet4.csv", skip=7, header=T) %>%
  mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_"),
                Sentrix_Row=substr(Sentrix_Position, 1, 3),
                Sentrix_Col=substr(Sentrix_Position, 4, 6)) %>%
  rename(shareid=Sample_Name) %>%
  select(shareid, sampleKey, Sentrix_ID, Sentrix_Row, Sentrix_Col)
# save("sampleData", file="../int/sampleData.RData")

### PHENOTYPE DATA
print("Phenotype variables...")

phenAll <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
  rename(sex=SEX, age=AGE8) %>%
  mutate(sex=as.character(factor(sex, levels=c(1,2), labels=c("M","F")))) %>%
  select(shareid, dbGaP_Subject_ID, sex, age, smoking_now, CVD_med, HT_med, T2D_med, SBP)
phen2 <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
  select(shareid, bmi)
labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
  select(-dbGaP_Subject_ID, -IDTYPE)

phenoData <- phenAll %>%
  left_join(phen2, by="shareid") %>%
  left_join(labValData, by="shareid")
source("helpers.R")  # Has Framingham Risk Score calculator function
phenoData$frs <- calc_FRS(phenoData)
# save("phenoData", file="../int/phenoData.RData")

### CVD DATA
print("Cardiovascular event data...")

soe <- read.delim("../data/events/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)

eventData <- soe %>%
  inner_join(labValData, by="shareid") %>%  # B/c has Exam 8 draw dates
  mutate(timeToEvent=DATE-drawdate) %>%  # Days between sample and event (NOTE: specific values here will change)
  rename(event=EVENT) %>%
  filter(event %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
  group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
  mutate(pastEvent=any(timeToEvent<0)) %>%  # For each individual, has he/she had an event before Exam 8?
  filter(timeToEvent>0) %>%  # Don't keep events that occurred prior to methylation collection
  slice(which.min(timeToEvent)) %>%  # Only keep first relevant event per person -- allow people with existing CVD for now
  ungroup() %>%
  select(shareid, event, drawdate, timeToEvent, pastEvent)
# save("eventData", file="../int/eventData.RData")  # This is just the events, so not all subjects included

survData <- left_join(sampleData, eventData, by="shareid") %>%  # Basis for including all relevant samples (not just those with events)
  mutate(event=!is.na(timeToEvent)) %>%  # event field becomes logical, indicating whether an event has occurred
  replace_na(list(timeToEvent=max(.$timeToEvent, na.rm=T))) %>%  # Replace missing times with max observed event time ## NOTE: this should be switched for an *actual* number
  select(shareid, event, timeToEvent, pastEvent)

### OUTPUTS FOR DOWNSTREAM USE
metaData <- Reduce(function(x,y) inner_join(x, y, by="shareid"), list(sampleData, survData, phenoData))
save("metaData", file="../int/metaData.RData")

sampleSheet <- select(metaData, one_of(names(sampleData), names(survData)), sex, age)
save("sampleSheet", file="../int/sampleSheet.RData")
