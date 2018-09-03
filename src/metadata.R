# Assemble metadata (demographics, CVD events, lab values, etc.) on subjects

suppressMessages(silent <- lapply(c("tidyverse"), library, character.only=T))

### SAMPLE DATA ###
print("Sample metadata...")

sampleData_fhs_c1 <- read_tsv("../data/fhs/meth/sample_attributes_c1.txt", skip=10)
sampleData_fhs_c2 <- read_tsv("../data/fhs/meth/sample_attributes_c2.txt", skip=10)
sampleData_fhs <- bind_rows(sampleData_fhs_c1, sampleData_fhs_c2) %>%
  mutate(subjID=as.integer(gsub("_[0-9]*", "", SAMPID)),  # These subject IDs (shareids) are in the 0-30,000 range
         slide=gsub("_[[:alnum:]]*", "", LABID),
         arrayPos=gsub("[0-9]*_", "", LABID),
         sentrixRow=substr(arrayPos, 1, 3),
         sentrixCol=substr(arrayPos, 4, 6),
         center=ifelse(grepl("^IL",Plate), "JHU", "UofMinn")) %>%
  rename(sampleKey=LABID,
         plate=Plate) %>%
  select(subjID, sampleKey, sentrixRow, sentrixCol, plate, center)
  
sampleData_whi_c1 <- read_tsv("../data/whi/meth/sample_attributes_c1.txt", skip=10)
sampleData_whi_c2 <- read_tsv("../data/whi/meth/sample_attributes_c2.txt", skip=10)
sampleData_whi <- bind_rows(sampleData_whi_c1, sampleData_whi_c2) %>%
  mutate(subjID=SUBJECT_ID,  # These subject IDs are in the 700,000-900,000 range
         sampleKey=paste(Methyl_Array, Array_pos, sep="_"),
         sentrixRow=substr(Array_pos, 1, 3),
         sentrixCol=substr(Array_pos, 4, 6)) %>%
  rename(plate=Plate,
         dnaPull=BA23_Pull_ID_DNA) %>%
  select(subjID, sampleKey, sentrixRow, sentrixCol, plate, dnaPull)

sampleData <- bind_rows(sampleData_fhs, sampleData_whi, .id="study") %>%
  mutate(study=case_when(study==1~"fhs", study==2~"whi"))

### PHENOTYPE DATA ###

## FHS
phenos_fhs_c1 <- read_tsv("../data/fhs/phen/phenos_fhs_c1.txt", skip=10)
phenos_fhs_c2 <- read_tsv("../data/fhs/phen/phenos_fhs_c2.txt", skip=10)
phenos_fhs <- bind_rows(phenos_fhs_c1, phenos_fhs_c2) %>%
  filter(shareid %in% sampleData_fhs$subjID) %>%
  rename(subjID=shareid, sex=SEX, age=AGE8, bmi=BMI8, 
         smk_now=CURRSMK8, cig_per_day=CPD8,  
         sbp=SBP8, glu=FASTING_BG8, chol=TC8, ldl=CALC_LDL8, hdl=HDL8, tg=TRIG8,
         ht_med=HRX8, lipid_med=LIPRX8, dm_med=DMRX8) %>%
  mutate(sex=ifelse(sex==1, "M", "F"),
         race="white",
         pack_years=) %>%
  mutate_at(c("ht_med","lipid_med","dm_med"), as.logical) %>%
  select(subjID, sex, age, race, bmi, smk_now, cig_per_day, 
         sbp, glu, chol, ldl, hdl, tg, ht_med, lipid_med, dm_med)

crp_fhs_c1 <- read_tsv("../data/fhs/phen/crp_c1.txt", skip=10)
crp_fhs_c2 <- read_tsv("../data/fhs/phen/crp_c2.txt", skip=10)
crp_fhs <- bind_rows(crp_fhs_c1, crp_fhs_c2) %>%
  rename(subjID=shareid, hscrp=crp)

questionnaire_fhs_c1 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c1.txt", skip=10)
questionnaire_fhs_c2 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c2.txt", skip=10)
questionnaire_fhs <- bind_rows(questionnaire_fhs_c1, questionnaire_fhs_c2) %>%
  rename(subjID=shareid, cig_start_age=H065) %>%
  select(subjID, cig_start_age)

phenoData_fhs <- phenos_fhs %>%
  left_join(crp_fhs, by="subjID") %>%
  left_join(questionnaire_fhs, by="subjID") %>%
  mutate(smk_py=ifelse(smk_now, cig_per_day/20*(age-cig_start_age), 0))  # Very rough -- assumes no one has quit

## WHI
basicData_whi_c1 <- read_tsv("../data/whi/phen/basic_whi_c1.txt", skip=10)
basicData_whi_c2 <- read_tsv("../data/whi/phen/basic_whi_c2.txt", skip=10)
basicData_whi <- bind_rows(basicData_whi_c1, basicData_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID) %>%
  rename(subjID=SUBJID, age=AGE, race=RACE) %>%
  mutate(sex="F",
         race=c("1"="amind","2"="asian","3"="black","4"="hispanic","5"="white","8"="other")[race]) %>%
  select(subjID, sex, age, race)

behaviorData_whi_c1 <- read_tsv("../data/whi/phen/behavior_c1.txt", skip=10)
behaviorData_whi_c2 <- read_tsv("../data/whi/phen/behavior_c2.txt", skip=10)
behaviorData_whi <- bind_rows(behaviorData_whi_c1, behaviorData_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID) %>%
  rename(subjID=SUBJID, smk_now=SMOKNOW, smk_py=PACKYRS) %>%
  select(subjID, smk_now, smk_py)

examData_whi_c1 <- read_tsv("../data/whi/phen/physical_whi_c1.txt", skip=10)
examData_whi_c2 <- read_tsv("../data/whi/phen/physical_whi_c2.txt", skip=10)
examData_whi <- bind_rows(examData_whi_c1, examData_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID) %>%
  rename(subjID=SUBJID, sbp=SYST, bmi=BMIX, daysSinceEnrollment=F80DAYS) %>%
  group_by(subjID) %>%
  slice(which.min(daysSinceEnrollment)) %>%  # Get only earliest exam per person
  ungroup() %>% 
  select(subjID, sbp, bmi, daysSinceEnrollment)

labData1_whi_c1 <- read_tsv("../data/whi/phen/labs_c1.txt", skip=10)
labData1_whi_c2 <- read_tsv("../data/whi/phen/labs_c2.txt", skip=10)
labData1_whi <- bind_rows(labData1_whi_c1, labData1_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID,
         COREVTYP==1) %>%  # Only care about blood samples from first year 
  rename(subjID=SUBJID, glu=COREGLUC, chol=CORETCHO, ldl=CORELDLC, hdl=COREHDLC, tg=CORETRI) %>%
  select(subjID, glu, chol, ldl, hdl, tg)

drawData2_whi_c1 <- read_tsv("../data/whi/phen/draws2_c1.txt", skip=10) %>%
  select(DRAWID, DRAWVTYP)
drawData2_whi_c2 <- read_tsv("../data/whi/phen/draws2_c2.txt", skip=10) %>%
  select(DRAWID, DRAWVTYP)
drawData2_whi <- bind_rows(drawData2_whi_c1, drawData2_whi_c2) %>%
  filter(DRAWVTYP==1) %>%  # Only care about blood samples from first year 
  select(DRAWID)
labData2_whi_c1 <- read_tsv("../data/whi/phen/labs2_c1.txt", skip=10)
labData2_whi_c2 <- read_tsv("../data/whi/phen/labs2_c2.txt", skip=10)
labData2_whi <- bind_rows(labData2_whi_c1, labData2_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID,
         SPECTYPE=="Serum",
         TESTABBR %in% c("LDLC","HDLC","TCHO","TRI","GLUC","CRP","INSU")) %>%
  inner_join(drawData2_whi, by="DRAWID") %>%
  group_by(SUBJID, TESTABBR) %>%
  summarise(TESTVAL=mean(TESTVAL, na.rm=T)) %>%
  spread(key=TESTABBR, value=TESTVAL) %>%
  rename(subjID=SUBJID, chol=TCHO, ldl=LDLC, hdl=HDLC, tg=TRI, glu=GLUC, hscrp=CRP, ins=INSU)

labData_whi <- bind_rows(labData1_whi, labData2_whi) %>%
  group_by(subjID) %>%
  summarise_all(mean, na.rm=T)  # Take mean when there are duplicate individuals from CORE and non-CORE

medsData_whi_c1 <- read_tsv("../data/whi/phen/medications_c1.txt", skip=10)
medsData_whi_c2 <- read_tsv("../data/whi/phen/medications_c2.txt", skip=10)
medsRef_whi <- read_tsv("../data/whi/phen/medication_classes.dat")
medsData_whi <- bind_rows(medsData_whi_c1, medsData_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID,
         F44VY==1) %>%
  inner_join(medsRef_whi, by="TCCODE") %>%
  mutate(ht_med=grepl("DIURETIC|CALCIUM BLOCKER|ACE INHIBITOR|ANGIOTENSIN II|BETA BLOCKER|BETA-BLOCKER|ALPHA 1|ALPHA-2|VASODILATOR|ALDOSTERONE", TCNAME),
         lipid_med=grepl("HMG COA REDUCTASE", TCNAME),
         dm_med=grepl("INSULIN|GLUCOSIDASE|BIGUANIDE|MEGLITINIDE|SULFONYLUREA|THIAZOLIDINEDIONES", TCNAME)) %>%
  rename(subjID=SUBJID) %>%
  group_by(subjID) %>%
  summarise(ht_med=any(ht_med),
            lipid_med=any(lipid_med),
            dm_med=any(dm_med)) %>%
  select(subjID, ht_med, lipid_med, dm_med)

phenoData_whi <- Reduce(function(x,y) left_join(x,y,by="subjID"), 
                        list(basicData_whi, behaviorData_whi, examData_whi, labData_whi, medsData_whi))

phenoData <- bind_rows(phenoData_fhs, phenoData_whi, .id="study") %>%
  mutate(study=case_when(study==1~"fhs", study==2~"whi"))


### CVD DATA ###
print("Cardiovascular event data...")

## FHS
exam_dates_fhs_c1 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c1.txt", skip=10, 
                              col_types=cols_only(shareid="i",date8="i",date9="i"))
exam_dates_fhs_c2 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c2.txt", skip=10,
                              col_types=cols_only(shareid="i",date8="i",date9="i"))
exam_dates_fhs <- bind_rows(exam_dates_fhs_c1, exam_dates_fhs_c2)

soe2015_fhs_c1 <- read_tsv("../data/fhs/phen/soe2015_c1.txt", skip=10, col_types=cols(shareid="i"))
soe2015_fhs_c2 <- read_tsv("../data/fhs/phen/soe2015_c2.txt", skip=10, col_types=cols(shareid="i"))
soe2015_fhs <- bind_rows(soe2015_fhs_c1, soe2015_fhs_c2)

soe2015_fhs_clean <- inner_join(soe2015_fhs, exam_dates_fhs, by="shareid") %>%
  filter(!is.na(date8),  # Had a visit during Exam 8 
         EVENT %in% c(1:29)) %>%  # CVD events (my definition) and all death
  mutate(eventType=case_when(EVENT %in% c(1:9,21:24) ~ "chd",  # MI variants, AP, CHD deaths
                             EVENT %in% c(10:19,25) ~ "stroke",  # CVA variants, ABI, TIA, embolism, hemorrhage
                             EVENT == 26 ~ "death_otherCVD",
                             EVENT %in% 27:29 ~ "death_nonCVD"),
         cvd=eventType %in% c("chd","stroke","death_otherCVD"),
         time=DATE-date8) %>%
  group_by(shareid) %>%
  summarise(pastEvent=any(cvd==T & time<=0),  # Note if subject had an event before Exam 8
            event=any(cvd==T & time>0),  # Future event if occurred after Exam 8
            timeToEvent=ifelse(any(cvd==T & time>0), min(time[cvd==T & time>0]), NA),  # Earliest post-Exam 8 event time
            eventType=ifelse(any(cvd==T & time>0), eventType[which.min(time[cvd==T & time>0])], as.character(NA)),
            death=any(EVENT %in% 21:29),  # Did the person die?
            timeToDeath=ifelse(any(death), time[EVENT %in% 21:29], NA),  # If they died, get time to death
            incCHD=any(eventType=="chd" & time>0),
            incStroke=any(eventType=="stroke" & time>0),
            timeToExam9=median(date9-date8))  # Carry through for censorship times

surv2014_fhs_c1 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c1.txt", skip=10)
surv2014_fhs_c2 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c2.txt", skip=10)
surv2014_fhs <- bind_rows(surv2014_fhs_c1, surv2014_fhs_c2)

outcomeData_fhs <- left_join(surv2014_fhs, soe2015_fhs_clean, by="shareid") %>%
  filter(shareid %in% sampleData_fhs$subjID) %>%
  replace_na(list(event=F, pastEvent=F, death=F)) %>%  # Add "false" for events/deaths after left join
  mutate(time=case_when(event==T ~ timeToEvent,  # Experienced an event after Exam 8 -> use that follow-up time
                        death==T ~ timeToDeath,  # Didn't experience an event but died -> censor at death
                        cvd==0 ~ cvddate,  # No event and no CVD in surv file -> censoring time from surv file
                        TRUE ~ as.integer(timeToExam9))) %>%  # No event and CVD in surv file -> use exam 9 date for censor time
  filter(!is.na(time)) %>%  # Remove those individuals with no Exam 9 date and thus no known censorship time
  rename(subjID=shareid) %>%
  select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)
## NOTE: THE APPROACH ABOVE LEAVES A NUMBER OF SUBJECTS (114) WHO HAD A PAST EVENT, NO FUTURE EVENT,
## AND NO AVAILABLE EXAM 9 DATE WITH MISSING "TIME" VALUES AND THEY ARE THUS EXCLUDED
## THIS IS POTENTIALLY REASONABLE HERE BECAUSE NONE REPRESENT CASES (THE MORE CRITICAL SAMPLES TO KEEP)

## WHI
whi_event_names <- c("CHD","CREVASC","STROKE")
whi_event_days <- paste0(whi_event_names, "DY")
outcomeData_whi_c1 <- read_tsv("../data/whi/phen/outcomes_ctos_c1.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
outcomeData_whi_c2 <- read_tsv("../data/whi/phen/outcomes_ctos_c2.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
outcomeData_whi <- bind_rows(outcomeData_whi_c1, outcomeData_whi_c2) %>%
  filter(SUBJID %in% sampleData_whi$subjID) %>%
  mutate(event=rowSums(.[,whi_event_names], na.rm=T)>0,
         time=do.call(pmin, c(.[,c(whi_event_days,"ENDFOLLOWDY")], na.rm=T)),  # Event time or censoring time
         eventType=case_when(time==CHDDY ~ "chd",
                             time==CREVASCDY ~ "crevasc",
                             time==STROKEDY ~ "stroke",
                             TRUE ~ as.character(NA)),
         incCHD=CHD==1,
         incStroke=STROKE==1,
         pastEvent=F) %>%  
  rename(subjID=SUBJID) %>%
  select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)

outcomeData <- bind_rows(outcomeData_fhs, outcomeData_whi, .id="study") %>%
  mutate(study=case_when(study==1~"fhs", study==2~"whi"))
  
# outcomeData_whi_c1_nonCaD <- read_tsv("../data/whi/phen/outcomes_cardio_c1.txt", skip=10)
# outcomeData_whi_c2_nonCaD <- read_tsv("../data/whi/phen/outcomes_cardio_c2.txt", skip=10)
# outcomeData_whi_nonCaD <- bind_rows(outcomeData_whi_c1_nonCaD, outcomeData_whi_c2_nonCaD) %>%
#   select(MI, CREVASC, MIDX, THRMAGNT)
# outcomeData_whi_c1_CaD <- read_tsv("../data/whi/phen/outcomes_cardio_cad_c1.txt", skip=10)
# outcomeData_whi_c2_CaD <- read_tsv("../data/whi/phen/outcomes_cardio_cad_c2.txt", skip=10)
# outcomeData_whi_CaD <- bind_rows(outcomeData_whi_c1_CaD, outcomeData_whi_c2_CaD) %>%
#   select(MI, CREVASC, MIDX, THRMAGNT) 
# outcomeData_whi <- bind_rows(outcomeData_whi_nonCaD, outcomeData_whi_CaD)


### OUTPUTS FOR DOWNSTREAM USE
metaData <- Reduce(function(x,y) inner_join(x, y),  # Default natural join by subjID and study 
                   list(sampleData, phenoData, outcomeData))
saveRDS(metaData, file="../int/metaData.rds")

sampleSheet <- select(metaData, one_of(names(sampleData)), sex, age, pastEvent, event)
saveRDS(sampleSheet, file="../int/sampleSheet.rds")



### PAST EXAM DATA FROM FHS ###
crp_ex2_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex2_c1.txt", skip=10, col_types=cols(shareid="i"))
crp_ex2_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex2_c2.txt", skip=10, col_types=cols(shareid="i"))
crp_ex2 <- bind_rows(crp_ex2_c1, crp_ex2_c2)
crp_ex6_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex6_c1.txt", skip=10, col_types=cols(shareid="i"))
crp_ex6_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex6_c2.txt", skip=10, col_types=cols(shareid="i"))
crp_ex6 <- bind_rows(crp_ex6_c1, crp_ex6_c2)
crp_ex7_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex7_c1.txt", skip=10, col_types=cols(shareid="i"))
crp_ex7_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex7_c2.txt", skip=10, col_types=cols(shareid="i"))
crp_ex7 <- bind_rows(crp_ex7_c1, crp_ex7_c2)
past_crp_data <- bind_rows(list(exam2=crp_ex2, exam6=crp_ex6, exam7=crp_ex7), .id="exam") %>%
  rename(hscrp=CRP,
         subjID=shareid) %>%
  select(subjID, hscrp, exam)

fram_lipid_varNames <- data.frame(exam=1:7,
                                  weight_lbs=c("A50","B15","C416","D401","E024","F007","G440"),
                                  height_inches=c("A51","B16","C417","D402","E025","F008","G441"),
                                  chol=c("A9","B352","C429","D448","E667","F726","G704"),
                                  hdl=c("A10","B355","C431","D449","E668","F725","G703"),
                                  ldl=c("A12","B357","C442",NA,NA,NA,NA),
                                  tg=c("A13","B358","C433","D451","E670","F727","G706"),
                                  glu=c("A31","B737","C434","D452","E671","F724","G705"),
                                  sbp1=c("A55","B24","C184","D192","E485","F476","G271"),
                                  sbp2=c("A57","B26","C286","D288","E581","F576","G354"))
read_past_exam <- function(row) {
  c1 <- read_tsv(paste0("../data/fhs/phen/past_exams/blood_ex", row["exam"], "_c1.txt"), 
           skip=10, col_types=cols(shareid="i", A24="i", B369="i"))
  c2 <- read_tsv(paste0("../data/fhs/phen/past_exams/blood_ex", row["exam"], "_c2.txt"), 
           skip=10, col_types=cols(shareid="i", A24="i", B369="i"))
  both <- bind_rows(c1, c2)
  cols <- c(shareid="shareid", na.omit(row)[-1])
  setNames(both[cols], names(cols))
}
past_exams_separate <- apply(fram_lipid_varNames, 1, read_past_exam)
names(past_exams_separate) <- paste0("exam", 1:7)
past_exam_data <- bind_rows(past_exams_separate, .id="exam") %>%
  mutate(bmi=weight_lbs/height_inches^2*703,
         sbp=rowMeans(.[,c("sbp1","sbp2")], na.rm=T),
         nonHDL=chol-hdl) %>%
  select(-one_of("weight_lbs","height_inches","sbp1","sbp2")) %>%
  dplyr::rename(subjID=shareid)

past_exam_data_withCrp <- full_join(past_exam_data, past_crp_data, by=c("exam","subjID"))
saveRDS(past_exam_data_withCrp, "../int/pastExamData.rds")