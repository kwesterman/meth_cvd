suppressMessages(silent <- lapply(c("tidyverse","parallel","survival","gtools"), library, character.only=T))

args <- commandArgs(trailingOnly=T)
job_idx <- as.numeric(args[1])  # Index in job array
total_jobs <- as.numeric(args[2])  # Total number of jobs in array
num_cores <- as.numeric(args[3])  # Get number of available cores as script argument

# # Load covariate data
# AgeSmokeData <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
#   dplyr::select(shareid, dbGaP_Subject_ID, SEX, AGE8, smoking_now)
# BMIData <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
#   dplyr::select(shareid, bmi)
# labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
#   dplyr::select(shareid, drawdate)
# 
# covData <- AgeSmokeData %>%
#   inner_join(BMIData, by="shareid") %>%
#   dplyr::rename(sex=SEX, age=AGE8)
# 
# # Load control probe PCs
# load("../int/CP_PCs.RData")
# CP_PCs <- rownames_to_column(data.frame(CP_PCs), var="sampleKey")
# 
# # Load cell counts
# cellCount_files <- grep("cellCounts_[0-9]+.RData", list.files("../int/"), value=T)
# cellCount_list <- lapply(paste0("../int/", cellCount_files), function(x) {load(x); cellCounts})
# cellCounts <- do.call(rbind, cellCount_list) %>%
#   data.frame() %>%
#   rownames_to_column(var="sampleKey")
# 
# # Load methylation and sample data
# load("../int/Mvals.RData")  # Loads Mvals object
# print("Methylation data loaded.")
# relevant_sites <- seq(from=floor((job_idx-1)/total_jobs*nrow(Mvals))+1,  # Sets "chunk" of sample indices to analyze
#                       to=floor(job_idx/total_jobs*nrow(Mvals)))
# Mvals <- Mvals[relevant_sites,]
# 
# sampleData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
#   mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_")) %>%
#   dplyr::rename(shareid=Sample_Name) %>%
#   inner_join(CP_PCs, by="sampleKey") %>%  # Add CPA control probe PCs
#   dplyr::select(shareid, sampleKey, one_of(paste0("PC",1:30,"_cp"))) %>%
#   inner_join(cellCounts, by="sampleKey")  # Add estimated cell counts
# 
# # Load event data
# soe2015 <- read.delim("../data/events/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)
# 
# # Clean event data and set up survival object
# soe2015.clean <- soe2015 %>%
#   inner_join(labValData, by="shareid") %>%
#   dplyr::mutate(timeToEvent=DATE-drawdate) %>%  # Days between sample and event (NOTE: specific values here will change)
#   dplyr::filter(timeToEvent>0) %>%  # Don't keep events that occurred prior to methylation collection
#   dplyr::filter(EVENT %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
#   group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
#   dplyr::slice(which.min(timeToEvent)) %>%  # Only keep first relevant event per person -- allow people with existing CVD for now
#   ungroup() %>%
#   dplyr::select(shareid, dbGaP_Subject_ID, EVENT, drawdate, timeToEvent)


# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")
relevant_sites <- seq(from=floor((job_idx-1)/total_jobs*nrow(Mvals))+1,  # Sets "chunk" of sample indices to analyze
                      to=floor(job_idx/total_jobs*nrow(Mvals)))
print(paste("Sample indices for this chunk:", min(relevant_sites), "to", max(relevant_sites)))
Mvals <- Mvals[relevant_sites,]

# Load all covariate data and prepare survival object
load("../int/phenoData.RData")
load("../int/sampleData.RData")
load("../int/eventData.RData")

survData <- sampleData %>%
  left_join(phenoData, by="shareid") %>%
  left_join(eventData, by="shareid") %>%
  dplyr::mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
         time=na.replace(timeToEvent, 1000)) %>%  ## NOTE: this 800 needs to be switched for an *actual* number
  dplyr::slice(match(colnames(Mvals), sampleKey))  # Re-order so rows match methylation columns
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor

# Sanity checks
print("Dimensions of survData:")
dim(survData)
print("Dimensions of M-value matrix:")
dim(Mvals)
stopifnot(ncol(Mvals)==nrow(survData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==survData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

# Formulas for models of interest
model_list <- list(
  model_unadj="survObj~meth",
  model_basic="survObj~meth+sex+age+smoking_now+bmi",
  model_wbc="survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran",
  model_10PC=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                       paste0("PC",1:10,"_cp",collapse="+")),
  model_20PC=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                       paste0("PC",1:20,"_cp",collapse="+")),
  model_30PC=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                       paste0("PC",1:30,"_cp",collapse="+"))
)

run_cox <- function(i, model_spec) {
  # Given an index (corresponding to a methylation site), bind that methylation data
  # to the covariate data and run Cox proportional hazards regression
  # Return vector of NAs if the regression fails
  modelData <- cbind(survData, meth=Mvals[i,])
  returnVal <- tryCatch({
    cox.fit <- coxph(as.formula(model_spec), data=modelData)
    c(CpG=rownames(Mvals)[i], summary(cox.fit)$coef['meth',c('coef','z','Pr(>|z|)')])
  }, error=function(e) {
    print(e)
    return(c(rownames(Mvals)[i], NA, NA, NA))
  }, warning=function(w) {
    print(w)
    return(c(rownames(Mvals)[i], NA, NA, NA))
  })
  returnVal
} 

p1 <- proc.time()
cl <- makeCluster(num_cores)
clusterExport(cl, varlist=c("survData","Mvals"))
print(proc.time()-p1)
clusterEvalQ(cl, library(survival))
resMats <- lapply(model_list, function(model) {
  res <- parLapply(cl, 1:nrow(Mvals), run_cox, model) 
  resMat <- do.call(rbind, res)
  resMat
})
# res <- parLapply(cl, 1:nrow(Mvals), run_cox, model_spec)
stopCluster(cl)
print(proc.time()-p1)

# resMat <- do.call(rbind, res)
save("resMats", file=paste0("../int/resMats_",job_idx,".RData"))
