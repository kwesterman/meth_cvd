suppressMessages(lapply(c("tidyverse","parallel","survival","gtools"), library, character.only=T))

args <- commandArgs(trailingOnly=T)
num_cores <- as.numeric(args[1])  # Get number of available cores as script argument

# Load covariate data
AgeSmokeData <- read.csv("../data/phenotypes/Off_Exam8_phen_cov_all.csv") %>%  # Has age and smoking status
  dplyr::select(shareid, dbGaP_Subject_ID, SEX, AGE8, smoking_now)
BMIData <- read.csv("../data/phenotypes/Off_Exam8_ex1_phen2.csv") %>%  # Has BMI
  dplyr::select(shareid, bmi)
labValData <- read.csv("../data/phenotypes/Off_ex8_lab_measures.csv") %>%  # Has blood draw date
  dplyr::select(shareid, drawdate)

covData <- AgeSmokeData %>%
  inner_join(BMIData, by="shareid") %>%
  dplyr::rename(sex=SEX, age=AGE8)

# Load control probe PCs
load("../int/CP_PCs.RData")
CP_PCs <- rownames_to_column(data.frame(CP_PCs), var="sampleKey")

# Load cell counts
cellCount_files <- grep("cellCounts_[0-9]+.RData", list.files("../int/"), value=T)
cellCount_list <- lapply(paste0("../int/", cellCount_files), function(x) {load(x); cellCounts})
cellCounts <- do.call(rbind, cellCount_list) %>%
  data.frame() %>%
  rownames_to_column(var="sampleKey")

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

sampleData <- read.csv("../data/sample_sheet4.csv", skip=7, header=T) %>%
  mutate(sampleKey=paste(Sentrix_ID, Sentrix_Position, sep="_")) %>%
  dplyr::rename(shareid=Sample_Name) %>%
  inner_join(CP_PCs, by="sampleKey") %>%  # Add CPA control probe PCs
  dplyr::select(shareid, sampleKey, one_of(paste0("PC",1:20,"_cp"))) %>%
  inner_join(cellCounts, by="sampleKey")  # Add estimated cell counts

# Load event data
soe2015 <- read.delim("../data/events/phs000007.v29.pht000309.v12.p10.c1.vr_soe_2015_a_0991s.HMB-IRB-MDS.txt", skip=10)

# Clean event data and set up survival object
soe2015.clean <- soe2015 %>%
  inner_join(labValData, by="shareid") %>%
  dplyr::mutate(timeToEvent=DATE-drawdate) %>%  # Days between sample and event (NOTE: specific values here will change)
  dplyr::filter(timeToEvent>0) %>%  # Don't keep events that occurred prior to methylation collection
  dplyr::filter(EVENT %in% 1:26) %>%  # Only events of interest (CVD, stroke, etc. -- see data dictionary)
  group_by(shareid) %>%  # To trim to earliest relevant event per person (what about people who have already had events?)
  dplyr::slice(which.min(timeToEvent)) %>%  # Only keep first relevant event per person -- allow people with existing CVD for now
  ungroup() %>%
  dplyr::select(shareid, dbGaP_Subject_ID, EVENT, drawdate, timeToEvent)

survData <- covData %>%
  left_join(soe2015.clean, by="shareid") %>%
  left_join(sampleData, by="shareid") %>%
  mutate(event=!is.na(timeToEvent),  # event == TRUE when a specific time of event exists, otherwise FALSE
         time=na.replace(timeToEvent, 1000)) %>%  ## NOTE: this 800 needs to be switched for an *actual* number
  slice(match(colnames(Mvals), sampleKey))  # Re-order so rows match methylation columns
survData$survObj <- Surv(time=survData$time, event=survData$event, type="right")  # all subjects without an event time get max. # days as time to censor
stopifnot(ncol(Mvals)==nrow(survData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==survData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

## Assess that everything is as expected
print("Dimensions of survData:")
dim(survData)
print("Dimensions of M-value matrix:")
dim(Mvals)
str(survData)

model_spec <- paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                     paste0("PC",1:20,"_cp",collapse="+"))

run_cox <- function(i) {
  modelData <- cbind(survData, meth=Mvals[i,])
  # survData$meth <- mVals[i,]  # Take the i'th CpG from the methylation matrix
  # cox.fit <- coxph(survObj~meth+CD8T+CD4T+NK+Bcell+Mono+Gran, data=survData)
  returnVal <- tryCatch({
    cox.fit <- coxph(as.formula(model_spec), data=modelData)
    c(CpG=rownames(Mvals)[i], summary(cox.fit)$coef['meth',c('coef','Pr(>|t|)')])
    # summary(cox.fit)$coefficients[2,'Pr(>|z|)']
    # cox.fit
  }, error=function(e) {
    return(c(rownames(Mvals)[i], NA, NA))
  }, warning=function(w) {
    return(c(rownames(Mvals)[i], NA, NA))
  })
  # print(returnVal)
  returnVal
} 
p1 <- proc.time()
# pVals <- lapply(1:nrow(Mvals[1:1000,]), run_cox)
cl <- makeCluster(num_cores)
clusterExport(cl, varlist=c("survData","Mvals","model_spec"))
print(proc.time()-p1)
clusterEvalQ(cl, library(survival))
res <- parLapply(cl, 1:nrow(Mvals), run_cox)
stopCluster(cl)
# pVals <- mclapply(1:nrow(Mvals[1:10000,]), run_cox, Mvals, mc.cores=num_cores)
print(proc.time()-p1)

resMat <- do.call(rbind, res)
write.csv(resMat, "../output/resMat.csv", row.names=F)
print(head(pVals))
print(sum(is.na(pVals[,3])))
# ggsave("../output/qqplot.jpg", qqp)
