# Assemble metadata (demographics, CVD events, lab values, etc.) on subjects

suppressMessages(silent <- lapply(c("tidyverse", "foreign"), library, character.only=T))

# SAMPLE DATA ------------------------------------------------------------------

print("Sample metadata...")

sample_data_fhs_c1 <- read_tsv("../data/fhs/meth/sample_attributes_c1.txt", skip=10)
sample_data_fhs_c2 <- read_tsv("../data/fhs/meth/sample_attributes_c2.txt", skip=10)
sample_data_fhs <- bind_rows(sample_data_fhs_c1, sample_data_fhs_c2) %>%
  mutate(subjID=gsub("_[0-9]*", "", SAMPID),  # These subject IDs (shareids) are in the 0-30,000 range
         slide=gsub("_[[:alnum:]]*", "", LABID),
         arrayPos=gsub("[0-9]*_", "", LABID),
         sentrixRow=substr(arrayPos, 1, 3),
         sentrixCol=substr(arrayPos, 4, 6),
         center=ifelse(grepl("^IL", Plate), "JHU", "UofMinn")) %>%
  rename(sampleKey=LABID,
         plate=Plate) %>%
  select(subjID, sampleKey, sentrixRow, sentrixCol, plate, center)

sample_data_whi_c1 <- read_tsv("../data/whi/meth/sample_attributes_c1.txt", skip=10)
sample_data_whi_c2 <- read_tsv("../data/whi/meth/sample_attributes_c2.txt", skip=10)
sample_data_whi <- bind_rows(sample_data_whi_c1, sample_data_whi_c2) %>%
  mutate(subjID=as.character(SUBJECT_ID),  # These subject IDs are in the 700,000-900,000 range
         sampleKey=paste(Methyl_Array, Array_pos, sep="_"),
         sentrixRow=substr(Array_pos, 1, 3),
         sentrixCol=substr(Array_pos, 4, 6)) %>%
  rename(plate=Plate,
         dnaPull=BA23_Pull_ID_DNA) %>%
  select(subjID, sampleKey, sentrixRow, sentrixCol, plate, dnaPull)

sample_data_lbc <- read_csv("../data/lbc/meth/target_QC_age_sex_date.csv") %>%
  mutate(subjID=as.character(ID_raw),
         sentrixRow=substr(pos, 1, 3),
         sentrixCol=substr(pos, 4, 6)) %>%
  rename(sampleKey=Basename,
         wave=WAVE) %>%
  select(subjID, sampleKey, sentrixRow, sentrixCol, plate, wave, cohort)
sample_data_lbc21 <- sample_data_lbc %>%
  filter(cohort == "LBC21") %>%
  select(-cohort)
sample_data_lbc36 <- sample_data_lbc %>%
  filter(cohort == "LBC36") %>%
  select(-cohort)

sample_data <- bind_rows(sample_data_fhs, sample_data_whi, 
                         sample_data_lbc21, sample_data_lbc36, .id="study") %>%
  mutate(study=case_when(study == 1 ~ "fhs", study == 2 ~ "whi", 
                         study == 3 ~ "lbc21", study == 4 ~ "lbc36"))

# PHENOTYPE DATA ---------------------------------------------------------------

## FHS
phenos_fhs_c1 <- read_tsv("../data/fhs/phen/phenos_fhs_c1.txt", skip=10)
phenos_fhs_c2 <- read_tsv("../data/fhs/phen/phenos_fhs_c2.txt", skip=10)
phenos_fhs <- bind_rows(phenos_fhs_c1, phenos_fhs_c2) %>%
  filter(shareid %in% sample_data_fhs$subjID) %>%
  rename(subjID=shareid, sex=SEX, age=AGE8, bmi=BMI8, 
         smk_now=CURRSMK8, cig_per_day=CPD8,  
         sbp=SBP8, glu=FASTING_BG8, chol=TC8, ldl=CALC_LDL8, hdl=HDL8, tg=TRIG8,
         ht_med=HRX8, lipid_med=LIPRX8, dm_med=DMRX8) %>%
  mutate(sex=ifelse(sex == 1, "M", "F"),
         race="white") %>%
  mutate_at(c("ht_med", "lipid_med", "dm_med"), as.logical) %>%
  select(subjID, sex, age, race, bmi, smk_now, cig_per_day, 
         sbp, glu, chol, ldl, hdl, tg, ht_med, lipid_med, dm_med)

crp_fhs_c1 <- read_tsv("../data/fhs/phen/crp_c1.txt", skip=10)
crp_fhs_c2 <- read_tsv("../data/fhs/phen/crp_c2.txt", skip=10)
crp_fhs <- bind_rows(crp_fhs_c1, crp_fhs_c2) %>%
  rename(subjID=shareid, hscrp=crp)

questionnaire_fhs_c1 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c1.txt", skip=10)
questionnaire_fhs_c2 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c2.txt", skip=10)
questionnaire_fhs <- bind_rows(questionnaire_fhs_c1, questionnaire_fhs_c2) %>%
  rename(subjID=shareid, smk_now=H062, cig_per_day=H063, cig_start_age=H065, 
         cig_stop_age=H066) %>%
  select(subjID, cig_start_age, cig_stop_age)

pheno_data_fhs <- phenos_fhs %>%
  left_join(crp_fhs, by="subjID") %>%
  left_join(questionnaire_fhs, by="subjID") %>%
  mutate(subjID=as.character(subjID),
         smk_py=case_when(smk_now == 0 & cig_stop_age == 0 ~ 0,
                          smk_now == 0 & cig_stop_age != 0 ~ 15 / 20 *  # Assumes ~ median cig/d for those who have quit
                            (cig_stop_age - cig_start_age),
                          smk_now == 1 ~ cig_per_day / 20 * (age - cig_start_age),
                          TRUE ~ as.numeric(NA)))
         # smk_py=ifelse(smk_now, cig_per_day / 20 * (age - cig_start_age), 0))  # Rough -- assumes no one has quit

## WHI
basic_data_whi_c1 <- read_tsv("../data/whi/phen/basic_whi_c1.txt", skip=10)
basic_data_whi_c2 <- read_tsv("../data/whi/phen/basic_whi_c2.txt", skip=10)
basic_data_whi <- bind_rows(basic_data_whi_c1, basic_data_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID) %>%
  rename(subjID=SUBJID, age=AGE, race=RACE) %>%
  mutate(sex="F",
         race=c("1"="amind", "2"="asian", "3"="black", "4"="hispanic", 
                "5"="white", "8"="other")[race]) %>%
  select(subjID, sex, age, race)

behavior_data_whi_c1 <- read_tsv("../data/whi/phen/behavior_c1.txt", skip=10)
behavior_data_whi_c2 <- read_tsv("../data/whi/phen/behavior_c2.txt", skip=10)
behavior_data_whi <- bind_rows(behavior_data_whi_c1, behavior_data_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID) %>%
  rename(subjID=SUBJID, smk_now=SMOKNOW, smk_py=PACKYRS) %>%
  select(subjID, smk_now, smk_py)

exam_data_whi_c1 <- read_tsv("../data/whi/phen/physical_whi_c1.txt", skip=10)
exam_data_whi_c2 <- read_tsv("../data/whi/phen/physical_whi_c2.txt", skip=10)
exam_data_whi <- bind_rows(exam_data_whi_c1, exam_data_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID) %>%
  rename(subjID=SUBJID, sbp=SYST, bmi=BMIX, daysSinceEnrollment=F80DAYS) %>%
  group_by(subjID) %>%
  slice(which.min(daysSinceEnrollment)) %>%  # Get only earliest exam per person
  ungroup() %>% 
  select(subjID, sbp, bmi, daysSinceEnrollment)

lab_data_1_whi_c1 <- read_tsv("../data/whi/phen/labs_c1.txt", skip=10)
lab_data_1_whi_c2 <- read_tsv("../data/whi/phen/labs_c2.txt", skip=10)
lab_data_1_whi <- bind_rows(lab_data_1_whi_c1, lab_data_1_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID,
         COREVTYP == 1) %>%  # Only care about blood samples from first year 
  rename(subjID=SUBJID, glu=COREGLUC, chol=CORETCHO, 
         ldl=CORELDLC, hdl=COREHDLC, tg=CORETRI) %>%
  select(subjID, glu, chol, ldl, hdl, tg)

drawData2_whi_c1 <- read_tsv("../data/whi/phen/draws2_c1.txt", skip=10) %>%
  select(DRAWID, DRAWVTYP)
drawData2_whi_c2 <- read_tsv("../data/whi/phen/draws2_c2.txt", skip=10) %>%
  select(DRAWID, DRAWVTYP)
drawData2_whi <- bind_rows(drawData2_whi_c1, drawData2_whi_c2) %>%
  filter(DRAWVTYP == 1) %>%  # Only care about blood samples from first year 
  select(DRAWID)
lab_data_2_whi_c1 <- read_tsv("../data/whi/phen/labs2_c1.txt", skip=10)
lab_data_2_whi_c2 <- read_tsv("../data/whi/phen/labs2_c2.txt", skip=10)
lab_data_2_whi <- bind_rows(lab_data_2_whi_c1, lab_data_2_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID,
         SPECTYPE == "Serum",
         TESTABBR %in% c("LDLC", "HDLC", "TCHO", "TRI", 
                         "GLUC", "CRP", "INSU")) %>%
  inner_join(drawData2_whi, by="DRAWID") %>%
  group_by(SUBJID, TESTABBR) %>%
  summarise(TESTVAL=mean(TESTVAL, na.rm=T)) %>%
  spread(key=TESTABBR, value=TESTVAL) %>%
  rename(subjID=SUBJID, chol=TCHO, ldl=LDLC, hdl=HDLC, tg=TRI, 
         glu=GLUC, hscrp=CRP, ins=INSU)

labData_whi <- bind_rows(lab_data_1_whi, lab_data_2_whi) %>%
  group_by(subjID) %>%
  summarise_all(mean, na.rm=T)  # Take mean when there are duplicate individuals from CORE and non-CORE

medsData_whi_c1 <- read_tsv("../data/whi/phen/medications_c1.txt", skip=10)
medsData_whi_c2 <- read_tsv("../data/whi/phen/medications_c2.txt", skip=10)
medsRef_whi <- read_tsv("../data/whi/phen/medication_classes.dat")
lipid_med_classes <- c(
  "ANTIHYPERLIPIDEMIC", "FIBRIC ACID DERIVATIVES",
  "INTESTINAL CHOLESTEROL ABSORPTION INHIBITORS", 
  "HMG COA REDUCTASE INHIBITORS", "HMG COA REDUCTASE INHIBITOR COMBINATIONS",
  "NICOTINIC ACID DERIVATIVES", "MISC. ANTIHYPERLIPIDEMICS", 
  "ANTIHYPERLIPIDEMIC COMBINATIONS", "FIBRIC ACID DERIVATIVE COMBINATIONS",
  "CALCIUM BLOCKER & HMG COA REDUCTASE INHIBITOR COMB"
)
diabetes_med_classes <- c(
  "ANTIDIABETIC", "INSULIN", "MIXED INSULIN", "BEEF INSULIN", "PORK INSULIN",
  "HUMAN INSULIN", "SULFONYLUREAS", "SULFONYLUREA COMBINATIONS",
  "ANTIDIABETIC - AMINO ACID DERIVATIVES", 
  "ANTIDIABETIC - D-PHENYLALANINE DERIVATIVES", "BIGUANIDES", 
  "MEGLITINIDE ANALOGUES", "DIABETIC OTHER", "DIABETIC OTHER - COMBINATIONS",
  "ALPHA-GLUCOSIDASE INHIBITORS", "INSULIN SENSITIZING AGENTS", 
  "THIAZOLIDINEDIONES", "ANTIDIABETIC COMBINATIONS",
  "SULFONYLUREA-BIGUANIDE COMBINATIONS",
  "THIAZOLIDINEDIONE-BIGUANIDE COMBINATIONS"
)
hypertension_med_classes <- c(
  "MINERALOCORTICOIDS", "VASOPRESSIN", "BETA BLOCKERS", 
  "BETA BLOCKERS NON-SELECTIVE", "BETA BLOCKERS CARDIO-SELECTIVE",
  "ALPHA-BETA BLOCKERS", "CALCIUM BLOCKERS", "ANTIHYPERTENSIVE",
  "ACE INHIBITORS", "ANGIOTENSIN II RECEPTOR ANTAGONIST",
  "ANTIADRENERGIC ANTIHYPERTENSIVES", "ANTIADRENERGICS - CENTRALLY ACTING",
  "ANTIADRENERGICS - PERIPHERALLY ACTING", "RESERPINE",
  "SELECTIVE ALDOSTERONE RECEPTOR ANTAGONISTS (SARAS)", "VASODILATORS",
  "FLUOROQUINOLONE VASODILATORS", "DOPAMINE D1 RECEPTOR AGONISTS",
  "ANTIHYPERTENSIVE - MAOIS", "MISC. ANTIHYPERTENSIVES",
  "ANTIHYPERTENSIVE COMBINATIONS", "RESERPINE COMBINATIONS",
  "ACE INHIBITORS & CALCIUM BLOCKERS", 
  "ACE INHIBITORS & THIAZIDE/THIAZIDE-LIKE",
  "BETA BLOCKER & DIURETIC COMBINATIONS", 
  "BETA BLOCKER & CALCIUM BLOCKER COMBINATIONS",
  "ANGIOTENSIN II RECEPTOR ANTAGONISTS & THIAZIDES",
  "ADRENOLYTICS-CENTRAL & THIAZIDE COMBINATIONS",
  "ADRENOLYTICS-PERIPHERAL & THIAZIDES", "ANTIHYPERTENSIVES-MAOIS & THIAZIDES",
  "ANTIHYPERTENSIVES-MISC & THIAZIDES", "VASODILATORS & THIAZIDES",
  "DIURETICS", "CARBONIC ANHYDRASE INHIBITORS", "LOOP DIURETICS",
  "MERCURIAL DIURETICS", "OSMOTIC DIURETICS", "POTASSIUM SPARING DIURETICS",
  "THIAZIDES AND THIAZIDE-LIKE DIURETICS", "MISCELLANEOUS DIURETICS",
  "COMBINATION DIURETICS", "DIURETICS & POTASSIUM", 
  "NON PRESCRIPTION DIURETICS", "PRESSORS", "PRESSOR COMBINATIONS",
  "PERIPHERAL VASODILATORS", "VASODILATOR COMBINATIONS",
  "MICROVASODILATORS", "PROSTAGLANDIN VASODILATORS",
  "VASOACTIVE NATRIURETIC PEPTIDES", "VASOCONSTRICTOR INHIBITORS",
  "CALCIUM BLOCKER & HMG COA REDUCTASE INHIBITOR COMB"
)

medsData_whi <- bind_rows(medsData_whi_c1, medsData_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID,
         F44VY == 1) %>%
  inner_join(medsRef_whi, by="TCCODE") %>%
  mutate(ht_med=TCNAME %in% hypertension_med_classes,
         lipid_med=TCNAME %in% lipid_med_classes,
         dm_med=TCNAME %in% diabetes_med_classes) %>%
  # mutate(ht_med=grepl(paste("DIURETIC", "CALCIUM BLOCKER", "ACE INHIBITOR", 
  #                           "ANGIOTENSIN II", "BETA BLOCKER", "BETA-BLOCKER", 
  #                           "ALPHA 1", "ALPHA-2", "VASODILATOR", "ALDOSTERONE", 
  #                           sep="|"),
  #                     TCNAME),
  #        lipid_med=grepl("HMG COA REDUCTASE", TCNAME),
  #        dm_med=grepl(paste("INSULIN", "GLUCOSIDASE", "BIGUANIDE", 
  #                           "MEGLITINIDE", "SULFONYLUREA", 
  #                           "THIAZOLIDINEDIONES", sep="|"),
  #                     TCNAME)) %>%
  rename(subjID=SUBJID) %>%
  group_by(subjID) %>%
  summarise(ht_med=any(ht_med),
            lipid_med=any(lipid_med),
            dm_med=any(dm_med)) %>%
  select(subjID, ht_med, lipid_med, dm_med)

pheno_data_whi <- Reduce(function(x,y) left_join(x,y,by="subjID"), 
                         list(basic_data_whi, behavior_data_whi, exam_data_whi, 
                              labData_whi, medsData_whi)) %>%
  mutate(subjID=as.character(subjID))

## LBC 1921
pheno_data_lbc21 <- read.spss("../data/lbc/phen/LBC1921_MethylationCardioVascularRiskMeditereanDiet_JO_TUFTS_02FEB2018.sav") %>%
  as.data.frame(stringsAsFactors=F) %>%
  mutate(subjID=trimws(studyno),
         sex=ifelse(gender == "male", "M", "F"),
         age_w1=mhtdays/365,
         age_w3=agedaywtcrf/365,
         age_w4=agedays_w4/365,
         age_w5=agedays_w5/365,
         age=floor(age_w1),
         race="white",
         smk_now=(smoker == "yes, current smoker") * 1,
         ht_med=as.logical(betabloc) | as.logical(aceinhib),
         lipid_med=as.logical(stat),
         ## diab_med?
         diabetes=grepl("yes", diab, ignore.case=T),
         chol=ifelse(chol == 9999, NA, chol * 38.67),  # mmol/L to mg/dL
         tg=ifelse(triglyc == 9999, NA, triglyc * 88.57),
         sbp=as.integer(as.character(sitsys))) %>%  # Weird values here (too high?)
  select(subjID, sex, age, race, bmi, smk_now, sbp, chol, tg, ht_med, lipid_med)

## LBC 1936
statin_strings <- c("statin", "zocor", "lipitor")
statin_regex <- paste(statin_strings, collapse="|")
other_lipid_med_strings <- c("ezetimibe")
other_lipid_med_regex <- paste(other_lipid_med_strings, collapse="|")
ht_med_strings <- c(
  "zide", "bumetanide", "amiloride", "emide",  # diuretics
  "diltiazem", "dipine",  # Ca channel blockers
  "pril",  # ACE inhibitors
  "artan",  # Angiotensin II inhibitors
  "lol",  # Beta blockers
  "zosin", "indoramin",  # Alpha blockers
  "lactone"  # Aldosterone receptor agonists
  # Alpha-2 receptor agonists
)
ht_med_regex <- paste(ht_med_strings, collapse="|")

pheno_data_lbc36 <- read.spss("../data/lbc/phen/LBC1936_MethylationCardioVascularRiskMeditereanDiet_JO_TUFTS_02FEB2018.sav") %>%
  as.data.frame(stringsAsFactors=F) %>%
  mutate(subjID=trimws(lbc36no),
         sex=ifelse(sex == "Male", "M", "F"),
         age=floor(agedays_w1/365),
         race="white",
         bmi=bmi_w1,
         smk_now=as.integer(smokcat_w1 == "current smoker"),
         sbp=rowMeans(.[c("sbp1sit_w1", "sbp2sit_w1", "sbp3sit_w1")], na.rm=T),
         chol=bld_choles_w1 * 38.67,
         hdl=bld_hdlchol_w1 * 38.67,
         tg=bld_triglyc_w1 * 88.57,
         hscrp=bld_crprot_w1,
         diabetes=diab_w1 == "Yes") %>%
  unite(all_drugs, matches("drug.+w1")) %>%
  mutate(statin=grepl(statin_regex, all_drugs, ignore.case=T),
         other_lipid_med=grepl(other_lipid_med_regex, all_drugs, ignore.case=T),
         lipid_med=statin | other_lipid_med,
         ht_med=grepl(ht_med_regex, all_drugs, ignore.case=T)) %>%
  select(subjID, sex, age, race, bmi, smk_now, sbp, 
         chol, hdl, tg, hscrp, diabetes,
         statin, lipid_med, ht_med, other_lipid_med)

pheno_data <- bind_rows(list(fhs=pheno_data_fhs, whi=pheno_data_whi, 
                             lbc21=pheno_data_lbc21, lbc36=pheno_data_lbc36), 
                        .id="study")


# CVD DATA ---------------------------------------------------------------------

print("Cardiovascular event data...")

## FHS
exam_dates_fhs_c1 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c1.txt", skip=10, 
                              col_types=cols_only(
                                shareid="i",date8="i",date9="i"))
exam_dates_fhs_c2 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c2.txt", skip=10,
                              col_types=cols_only(
                                shareid="i",date8="i",date9="i"))
exam_dates_fhs <- bind_rows(exam_dates_fhs_c1, exam_dates_fhs_c2)

soe2015_fhs_c1 <- read_tsv("../data/fhs/phen/soe2015_c1.txt", skip=10,
                           col_types=cols(shareid="i"))
soe2015_fhs_c2 <- read_tsv("../data/fhs/phen/soe2015_c2.txt", skip=10, 
                           col_types=cols(shareid="i"))
soe2015_fhs <- bind_rows(soe2015_fhs_c1, soe2015_fhs_c2)

soe2015_fhs_clean <- inner_join(soe2015_fhs, exam_dates_fhs, by="shareid") %>%
  filter(!is.na(date8),  # Had a visit during Exam 8 
         EVENT %in% c(1:29)) %>%  # CVD events (my definition) and all death
  mutate(eventType=case_when(EVENT %in% c(1:9,21:24) ~ "chd",  # MI variants, AP, CHD deaths
                             EVENT %in% c(10:19,25) ~ "stroke",  # CVA variants, ABI, TIA, embolism, hemorrhage
                             EVENT == 26 ~ "death_otherCVD",
                             EVENT %in% 27:29 ~ "death_nonCVD"),
         cvd=eventType %in% c("chd", "stroke", "death_otherCVD"),
         time=DATE - date8) %>%
  group_by(shareid) %>%
  summarise(pastEvent=any(cvd == T & time <= 0),  # Note if subject had an event before Exam 8
            event=any(cvd == T & time > 0),  # Future event if occurred after Exam 8
            timeToEvent=ifelse(any(cvd == T & time > 0), 
                               min(time[cvd == T & time > 0]), NA),  # Earliest post-Exam 8 event time
            eventType=ifelse(any(cvd == T & time > 0), 
                             eventType[which.min(time[cvd == T & time > 0])], 
                             as.character(NA)),
            death=any(EVENT %in% 21:29),  # Did the person die?
            timeToDeath=ifelse(any(death), time[EVENT %in% 21:29], NA),  # If they died, get time to death
            incCHD=any(eventType == "chd" & time > 0),
            incStroke=any(eventType == "stroke" & time > 0),
            timeToExam9=median(date9 - date8))  # Carry through for censorship times

surv2014_fhs_c1 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c1.txt", skip=10)
surv2014_fhs_c2 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c2.txt", skip=10)
surv2014_fhs <- bind_rows(surv2014_fhs_c1, surv2014_fhs_c2)

outcome_data_fhs <- left_join(surv2014_fhs, soe2015_fhs_clean, by="shareid") %>%
  inner_join(exam_dates_fhs, by="shareid") %>%
  filter(shareid %in% sample_data_fhs$subjID) %>%
  replace_na(list(event=F, pastEvent=F, death=F)) %>%  # Add "false" for events/deaths after left join
  mutate(time=case_when(event == T ~ timeToEvent,  # Experienced an event after Exam 8 -> use that follow-up time
                        death == T ~ timeToDeath,  # Didn't experience an event but died -> censor at death
                        cvd == 0 ~ cvddate - date8,  # No event and no CVD in surv file -> censoring time from surv file
                        TRUE ~ timeToExam9),  # No event and CVD in surv file -> use exam 9 date for censor time
         time=as.integer(time)) %>%  
  filter(!is.na(time)) %>%  # Remove those individuals with no Exam 9 date and thus no known censorship time
  mutate(subjID=as.character(shareid)) %>%
  select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)
## NOTE: THE APPROACH ABOVE LEAVES A NUMBER OF SUBJECTS (114) WHO HAD A PAST EVENT, NO FUTURE EVENT,
## AND NO AVAILABLE EXAM 9 DATE WITH MISSING "TIME" VALUES AND THEY ARE THUS EXCLUDED
## THIS IS POTENTIALLY REASONABLE HERE BECAUSE NONE REPRESENT CASES (THE MORE CRITICAL SAMPLES TO KEEP)

## WHI
whi_event_names <- c("CHD", "CREVASC", "STROKE")
whi_event_days <- paste0(whi_event_names, "DY")
outcome_data_whi_c1 <- read_tsv("../data/whi/phen/outcomes_ctos_c1.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
outcome_data_whi_c2 <- read_tsv("../data/whi/phen/outcomes_ctos_c2.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
outcome_data_whi <- bind_rows(outcome_data_whi_c1, outcome_data_whi_c2) %>%
  filter(SUBJID %in% sample_data_whi$subjID) %>%
  mutate(event=rowSums(.[,whi_event_names], na.rm=T)>0,
         time=do.call(pmin, c(.[,c(whi_event_days,"ENDFOLLOWDY")], na.rm=T)),  # Event time or censoring time
         eventType=case_when(time == CHDDY ~ "chd",
                             time == CREVASCDY ~ "crevasc",
                             time == STROKEDY ~ "stroke",
                             TRUE ~ as.character(NA)),
         incCHD=(CHD == 1),
         incStroke=(STROKE == 1),
         pastEvent=F) %>%  
  mutate(subjID=as.character(SUBJID)) %>%
  select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)

## LBC 1921
outcome_data_lbc21 <- read.spss("../data/lbc/phen/LBC1921_MethylationCardioVascularRiskMeditereanDiet_JO_TUFTS_02FEB2018.sav") %>%
  as.data.frame(stringsAsFactors=F) %>%
  mutate(subjID=trimws(studyno),
         pastEvent=(grepl("yes", vaschist, ignore.case=T) | cerbvasc == "yes"),
         event=(cvhist83 == "yes" | crvhist83 == "yes" | cvhist87 == "yes" | 
                  crvhist87 == "yes" | cvhist_w4 == "yes" | 
                  crvhist_w4 == "yes" | cvdhist_w5 == "Yes" | 
                  stroke_w5 == "Yes"),
         event=ifelse(is.na(event), F, event),
         time=case_when(cvhist83 == "yes" | crvhist83 == "yes" ~ 4 * 365,
                        cvhist87 == "yes" | crvhist87 == "yes" ~ 
                          agedaywtcrf-mhtdays,
                        cvhist_w4 == "yes" | crvhist_w4 == "yes" ~ 
                          agedays_w4-mhtdays,
                        cvdhist_w5 == "Yes" | stroke_w5 == "Yes" ~ 
                          agedays_w5-mhtdays,
                        !is.na(agedays_death) ~ agedays_death - 79,
                        TRUE ~ 14 * 365)) %>%
  select(subjID, pastEvent, event, time)

## LBC36
outcome_data_lbc36 <- read.spss("../data/lbc/phen/LBC1936_MethylationCardioVascularRiskMeditereanDiet_JO_TUFTS_02FEB2018.sav") %>%
  as.data.frame(stringsAsFactors=F) %>%
  mutate(subjID=trimws(lbc36no),
         pastEvent=cvdhist_w1 == "Yes" | stroke_w1 == "Yes",
         event=(cvdhist_w2 == "Yes" | stroke_w2 == "Yes" | 
                  cvdhist_w3 == "Yes" | stroke_w3 == "Yes" | 
                  cvdhist_w4 == "Yes" | stroke_w4 == "Yes"),
         incStroke=(stroke_w2 == "Yes" | 
                      (cvdhist_w3 == "No" & stroke_w3 == "Yes") | 
                      (cvdhist_w3 == "No" & cvdhist_w4 == "No" & stroke_w4 == "Yes")),
         event=ifelse(is.na(event), F, event),
         time=case_when(cvdhist_w2 == "Yes" | stroke_w2 == "Yes" ~ 4 * 365,
                        cvdhist_w3 == "Yes" | stroke_w3 == "Yes" ~ 8 * 365,
                        cvdhist_w4 == "Yes" | stroke_w4 == "Yes" ~ 11 * 365,
                        !is.na(agedays_death) ~ agedays_death - agedays_w1,
                        TRUE ~ 14 * 365)) %>%
  select(subjID, pastEvent, event, time, incStroke)

outcome_data <- bind_rows(outcome_data_fhs, outcome_data_whi,
                          outcome_data_lbc21, outcome_data_lbc36,
                          .id="study") %>%
  mutate(study=case_when(study == 1 ~ "fhs", study == 2 ~ "whi", 
                         study == 3 ~ "lbc21", study == 4 ~ "lbc36"))


# OUTPUTS FOR DOWNSTREAM USE ---------------------------------------------------

metaData <- Reduce(function(x,y) inner_join(x, y),  # Default natural join by subjID and study 
                   list(sample_data, pheno_data, outcome_data))
saveRDS(metaData, file="../int/metaData.rds")
write_csv(metaData, "../int/metaData.csv")

sampleSheet <- inner_join(sample_data, select(metaData, subjID, sex, age, pastEvent, event), by="subjID") %>%
  mutate(study=ifelse(grepl("lbc", study), "lbc", study))
saveRDS(sampleSheet, file="../int/sampleSheet.rds")


# PAST EXAM DATA FROM FHS ------------------------------------------------------

crp_ex2_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex2_c1.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex2_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex2_c2.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex2 <- bind_rows(crp_ex2_c1, crp_ex2_c2)
crp_ex6_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex6_c1.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex6_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex6_c2.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex6 <- bind_rows(crp_ex6_c1, crp_ex6_c2)
crp_ex7_c1 <- read_tsv("../data/fhs/phen/past_exams/crp_ex7_c1.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex7_c2 <- read_tsv("../data/fhs/phen/past_exams/crp_ex7_c2.txt", 
                       skip=10, col_types=cols(shareid="i"))
crp_ex7 <- bind_rows(crp_ex7_c1, crp_ex7_c2)
past_crp_data <- bind_rows(list(exam2=crp_ex2, exam6=crp_ex6, exam7=crp_ex7), 
                           .id="exam") %>%
  rename(hscrp=CRP,
         subjID=shareid) %>%
  select(subjID, hscrp, exam)

fram_lipid_varNames <- data.frame(
  exam=1:7,
  weight_lbs=c("A50", "B15", "C416", "D401", "E024", "F007", "G440"),
  height_inches=c("A51", "B16", "C417", "D402", "E025", "F008", "G441"),
  chol=c("A9", "B352", "C429", "D448", "E667", "F726", "G704"),
  hdl=c("A10", "B355", "C431", "D449", "E668", "F725", "G703"),
  ldl=c("A12", "B357", "C442", NA, NA, NA, NA),
  tg=c("A13", "B358", "C433", "D451", "E670", "F727", "G706"),
  glu=c("A31", "B737", "C434", "D452", "E671", "F724", "G705"),
  sbp1=c("A55", "B24", "C184", "D192", "E485", "F476", "G271"),
  sbp2=c("A57", "B26", "C286", "D288", "E581", "F576", "G354"))
read_past_exam <- function(row) {
  c1 <- read_tsv(paste0("../data/fhs/phen/past_exams/blood_ex", 
                        row["exam"], "_c1.txt"), 
                 skip=10, col_types=cols(shareid="i", A24="i", B369="i"))
  c2 <- read_tsv(paste0("../data/fhs/phen/past_exams/blood_ex", 
                        row["exam"], "_c2.txt"), 
                 skip=10, col_types=cols(shareid="i", A24="i", B369="i"))
  both <- bind_rows(c1, c2)
  cols <- c(shareid="shareid", na.omit(row)[-1])
  setNames(both[cols], names(cols))
}
past_exams_separate <- apply(fram_lipid_varNames, 1, read_past_exam)
names(past_exams_separate) <- paste0("exam", 1:7)
past_exam_data <- bind_rows(past_exams_separate, .id="exam") %>%
  mutate(bmi=weight_lbs / height_inches ^ 2 * 703,
         sbp=rowMeans(.[, c("sbp1", "sbp2")], na.rm=T),
         nonHDL=chol - hdl) %>%
  select(-one_of("weight_lbs", "height_inches", "sbp1", "sbp2")) %>%
  dplyr::rename(subjID=shareid)

past_exam_data_with_crp <- full_join(past_exam_data, past_crp_data, 
                                     by=c("exam", "subjID")) %>%
  mutate(subjID=as.character(subjID))
saveRDS(past_exam_data_with_crp, "../int/past_exam_data.rds")